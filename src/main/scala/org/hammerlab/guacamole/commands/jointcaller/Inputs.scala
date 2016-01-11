package org.hammerlab.guacamole.commands.jointcaller

import org.apache.http.client.utils.URLEncodedUtils

import scala.collection.JavaConversions
import scala.collection.mutable.ArrayBuffer

object Inputs {

  object TissueType extends Enumeration {
    val Normal = Value("normal")
    val Tumor = Value("tumor")
  }

  object Analyte extends Enumeration {
    val DNA = Value("dna")
    val RNA = Value("rna")
  }

  case class Input(name: String, path: String, tissueType: TissueType.Value, analyte: Analyte.Value) {
    override def toString: String = {
      "<Input '%s' of %s %s at %s >".format(name, tissueType, analyte, path)
    }

  }
  object Input {
    def apply(url: String, defaults: Option[Input] = None): Input = {
      val parsed = new java.net.URI(url)
      val urlWithoutFragment = new java.net.URI(
        parsed.getScheme,
        parsed.getUserInfo,
        parsed.getHost,
        parsed.getPort,
        parsed.getPath,
        parsed.getQuery,
        "").toString.stripSuffix("#") // set fragment to the empty string

      val keyValues = URLEncodedUtils.parse(parsed.getFragment, org.apache.http.Consts.UTF_8)
      var tissueType: Option[TissueType.Value] = defaults.map(_.tissueType)
      var analyte: Option[Analyte.Value] = defaults.map(_.analyte)
      var name = defaults.map(_.name).filter(_.nonEmpty).getOrElse(
        urlWithoutFragment.split('/').last.stripSuffix(".bam"))
      JavaConversions.iterableAsScalaIterable(keyValues).foreach(pair => {
        val value = pair.getValue.toLowerCase
        pair.getName.toLowerCase match {
          case "tissue_type" => tissueType = Some(TissueType.withName(value))
          case "analyte"     => analyte = Some(Analyte.withName(value))
          case "name"        => name = value
          case other => {
            throw new IllegalArgumentException(
              "Unsupported input property: %s in %s. Valid properties are: tissue_type, analyte, name".format(
                other, url))
          }
        }
      })
      if (tissueType.isEmpty) {
        throw new IllegalArgumentException("No tissue_type specified for %s".format(url))
      }
      if (analyte.isEmpty) {
        throw new IllegalArgumentException("No analyte specified for %s".format(url))
      }
      new Input(name, urlWithoutFragment, tissueType.get, analyte.get)
    }
  }

  def parseMultiple(urls: Seq[String]): Seq[Input] = {
    val inputs = ArrayBuffer.newBuilder[Input]
    var default = Input("", "", TissueType.Normal, Analyte.DNA)
    urls.map(url => {
      val result = Input(url, Some(default))
      default = Input("", "", TissueType.Tumor, Analyte.DNA)
      result
    })
  }

  case class GroupedInputs(normalDNA: Seq[Int], normalRNA: Seq[Int], tumorDNA: Seq[Int], tumorRNA: Seq[Int]) {
    def relative(indices: Seq[Int]): GroupedInputs = {
      val oldToNew = indices.zipWithIndex.toMap
      def transform(seq: Seq[Int]): Seq[Int] = seq.flatMap(oldToNew.get(_))
      GroupedInputs(transform(normalDNA), transform(normalRNA), transform(tumorDNA), transform(tumorRNA))
    }
  }
  object GroupedInputs {
    def apply(inputs: Seq[Input]): GroupedInputs = {
      def indices(function: Input => Boolean) = {
        inputs.zipWithIndex.filter(pair => function(pair._1)).map(_._2)
      }
      GroupedInputs(
        normalDNA = indices(input => input.tissueType == TissueType.Normal && input.analyte == Analyte.DNA),
        normalRNA = indices(input => input.tissueType == TissueType.Normal && input.analyte == Analyte.DNA),
        tumorDNA = indices(input => input.tissueType == TissueType.Tumor && input.analyte == Analyte.DNA),
        tumorRNA = indices(input => input.tissueType == TissueType.Tumor && input.analyte == Analyte.RNA))
    }
  }

}
