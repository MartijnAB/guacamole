/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hammerlab.guacamole.commands

import java.io.InputStreamReader
import java.text.{ DecimalFormatSymbols, DecimalFormat }
import java.util.Locale
import javax.script._

import com.google.common.io.CharStreams
import htsjdk.samtools.util.Locus
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import org.apache.commons.io.IOUtils
import org.apache.commons.lang.StringEscapeUtils
import org.apache.hadoop.conf.Configuration
import org.apache.hadoop.fs.{ Path, FileSystem }
import org.apache.http.client.utils.URLEncodedUtils
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.hammerlab.guacamole.Common.Arguments.NoSequenceDictionary
import org.hammerlab.guacamole.likelihood.Likelihood
import org.hammerlab.guacamole.pileup._
import org.hammerlab.guacamole.reads._
import org.hammerlab.guacamole._
import org.hammerlab.guacamole.variants.{Genotype, AlleleEvidence, CalledAllele}
import org.kohsuke.args4j.spi.StringArrayOptionHandler
import org.kohsuke.args4j.{ Option => Args4jOption, Argument }

import scala.collection.JavaConversions
import scala.collection.mutable.ArrayBuffer

object SomaticJoint {
  protected class Arguments extends DistributedUtil.Arguments with NoSequenceDictionary {

    @Argument(required = true, multiValued = true,
      usage = "FILE1 FILE2 FILE3")
    var inputs: Array[String] = Array.empty

    @Args4jOption(name = "--out-small-somatic-variants", usage = "Output path. Default: stdout")
    var outSmallSomaticVariants: String = ""

    @Args4jOption(name = "-q", usage = "Quiet: less stdout")
    var quiet: Boolean = false
  }

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
      var name = defaults.map(_.name).filter(_.nonEmpty).getOrElse(urlWithoutFragment)
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

    def parseMultiple(urls: Seq[String]): Seq[Input] = {
      val inputs = ArrayBuffer.newBuilder[Input]
      var default = Input("", "", TissueType.Normal, Analyte.DNA)
      urls.map(url => {
        val result = Input(url, Some(default))
        default = Input("", "", TissueType.Tumor, Analyte.DNA)
        result
      })
    }
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

  object VariantKind extends Enumeration {
    val DeNovo = Value("denovo")
    val LossOfHeterozygosity = Value("LossOfHeterozygosity")

  }

  object PhaseRelationship extends Enumeration {
    val equal = Value("equal")
    val superclone = Value("superclone")
    val subclone = Value("subclone")
  }

  case class Locus(contig: String, start: Int, length: Int) {}

  case class Phasing(
                      relationship: PhaseRelationship.Value,
                      locus1: Locus,
                      locus2: Locus,
                      allele1: String,
                      allele2: String,
                      supportingReads: Int) {}

  trait SmallVariant {
    val locus: Locus
    val ref: String
    val alt: String
  }

  case class GermlineSmallVariant(locus: Locus, ref: String, alt: String) extends SmallVariant {}

  case class SomaticSmallVariant(kind: VariantKind.Value, locus: Locus, ref: String, germline: Seq[String], alt: String) {}


  object Caller extends SparkCommand[Arguments] {
    override val name = "somatic-joint"
    override val description = "somatic caller for any number of samples from the same patient"

    override def run(args: Arguments, sc: SparkContext): Unit = {
      val inputs = Input.parseMultiple(args.inputs)

      if (!args.quiet) {
        println("Running on %d inputs:".format(inputs.length))
        inputs.foreach(input => println(input))
      }

      val readSets = inputs.zipWithIndex.map({
        case (input, index) => ReadSet(
          sc,
          input.path, false, Read.InputFilters.empty, token = index, contigLengthsFromDictionary = !args.noSequenceDictionary)
      })

      assert(readSets.forall(_.sequenceDictionary == readSets(0).sequenceDictionary),
        "Samples have different sequence dictionaries: %s."
          .format(readSets.map(_.sequenceDictionary.toString).mkString("\n")))

      val loci = Common.loci(args, readSets(0))
      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(args, loci, readSets.map(_.mappedReads): _*)

      val groupedInputs = GroupedInputs(inputs)

      // TODO: we currently are re-shuffling on every pileupFlatMap call.

      // Call germline small variants.

      val germlineCalls = DistributedUtil.pileupFlatMapMultipleRDDs(
        groupedInputs.normalDNA.map(readSets(_).mappedReads),
        lociPartitions,
        true,  // skip empty
        pileups => findGermlineVariants(
          groupedInputs.normalDNA.map(inputs(_)),
          groupedInputs.relative(groupedInputs.normalDNA),
          pileups))
    }
  }
  def combinedPileup(pileups: Seq[Pileup]) = {
    val elements = pileups.flatMap(_.elements)
    Pileup(pileups(0).referenceName, pileups(0).locus, pileups(0).referenceBase, elements)
  }
  def findGermlineVariants(inputs: Seq[Input],
                           groupedInputs: GroupedInputs,
                           pileups: Seq[Pileup] ): Iterator[variants.Genotype] = {

    val pileup = combinedPileup(groupedInputs.normalDNA.map(pileups(_)))

    val alleles = (Bases.pileup.referenceBase )
      distinctAlleles.filter(allele => allele.altBases.forall((Bases.isStandardBase _)))
    val allelesWithRef = if (alleles.contains(Allele()))
    val genotypes = for {
      i <- 0 until alleles.size
      j <- i until alleles.size
    } yield Genotype(alleles(i), alleles(j))
    val likelihoods = likelihoodsOfGenotypes(pileup.elements, genotypes, probabilityCorrect, prior, logSpace, normalize)
    genotypes.zip(likelihoods)



    val genotypeLikelihoods = Likelihood.likelihoodsOfAllPossibleGenotypesFromPileup(
      combinedNormalPileup,
      probabilityCorrect = Likelihood.probabilityCorrectIncludingAlignment,
      logSpace = false,
      normalize = false)

    val (genotype, probability) = genotypeLikelihoods.maxBy(_._2)
    if (genotype.hasVariantAllele) Iterator(genotype) else Iterator.empty
  }
}

