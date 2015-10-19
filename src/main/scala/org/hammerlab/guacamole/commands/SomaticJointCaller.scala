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
import java.util
import java.util.Locale
import javax.script._

import com.google.common.io.CharStreams
import htsjdk.samtools.util.Locus
import htsjdk.variant.variantcontext.{GenotypeBuilder, Allele, VariantContextBuilder, VariantContext}
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf._
import org.apache.commons.io.IOUtils
import org.apache.commons.lang.StringEscapeUtils
import org.apache.hadoop.conf.Configuration
import org.apache.hadoop.fs.{ Path, FileSystem }
import org.apache.http.client.utils.URLEncodedUtils
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.util.PhredUtils
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
  class Arguments extends DistributedUtil.Arguments with NoSequenceDictionary {

    @Argument(required = true, multiValued = true,
      usage = "FILE1 FILE2 FILE3")
    var inputs: Array[String] = Array.empty

    @Args4jOption(name = "--out-small-germline-variants", usage = "Output path. Default: stdout")
    var outSmallGermlineVariants: String = ""

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

  case class Phasing(
                      relationship: PhaseRelationship.Value,
                      locus1: Locus,
                      locus2: Locus,
                      allele1: String,
                      allele2: String,
                      supportingReads: Int) {}


  case class GermlineSmallVariant(
                                   contig: String,
                                   position: Long,
                                   ref: String,
                                   alts: Seq[String],
                                   likelihoods: Map[Seq[String], Double]) {}

  //case class SomaticSmallVariant(kind: VariantKind.Value, locus: Locus, ref: String, germline: Seq[String], alt: String) {}


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
          pileups)).collect

      Common.progress("Called %,d germline variants.".format(germlineCalls.length))

      val sampleName = readSets(0).reads.take(1)(0).sampleName

      if (args.outSmallGermlineVariants.nonEmpty) {
        Common.progress("Writing germline variants.")
        val codec = new VCFCodec()
        val writer = new VariantContextWriterBuilder()
          .setOutputFile(args.outSmallGermlineVariants)
          .setReferenceDictionary(readSets(0).sequenceDictionary.get.toSAMSequenceDictionary)
          .build
        val headerLines = new util.HashSet[VCFHeaderLine]()
        headerLines.add(new VCFHeaderLine("GT", "Genotype"))
        // headerLines.add(new VCFHeaderLine("FORMAT", "Format"))
        // headerLines.add(
        //   new VCFFormatHeaderLine("FT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Genotype filters."))
        headerLines.add(
          new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype."))

        // headerLines.add(new VCFInfoHeaderLine())
        val header = new VCFHeader(headerLines, JavaConversions.seqAsJavaList(Seq(sampleName)))
        header.setSequenceDictionary(readSets(0).sequenceDictionary.get.toSAMSequenceDictionary)
        header.setWriteCommandLine(true)
        header.setWriteEngineHeaders(true)
        writer.writeHeader(header)

        germlineCalls.foreach(call => {
          assert(call.alts.nonEmpty)

          val genotypeAlleles = if (call.alts.length == 1) {
            Seq(Allele.create(call.ref, true)) ++ call.alts.map(Allele.create(_, false))
          } else {
            call.alts.map(Allele.create(_, false))
          }
          val allAlleles = Seq(Allele.create(call.ref, true)) ++ call.alts.distinct.map(Allele.create(_, false))

          val genotype = new GenotypeBuilder(sampleName)
            .alleles(JavaConversions.seqAsJavaList(genotypeAlleles))
            .make
          val context = new VariantContextBuilder()
            .chr(call.contig)
            .start(call.position + 1)  // one based based inclusive
            .stop(call.position + 1 + math.max(call.ref.length - 1, 0))
            .genotypes(genotype)
            .alleles(JavaConversions.seqAsJavaList(allAlleles))
            .make
          writer.add(context)
        })
        writer.close()
      }
      Common.progress("Wrote %,d calls to %s".format(germlineCalls.length, args.outSmallGermlineVariants))
    }
  }

  def combinedPileup(pileups: Seq[Pileup]) = {
    val elements = pileups.flatMap(_.elements)
    Pileup(pileups(0).referenceName, pileups(0).locus, pileups(0).referenceBase, elements)
  }
  def logLikelihoodPileup(elements: Iterable[PileupElement], mixture: Map[String, Double]): Double = {
    def logLikelihoodPileupElement(element: PileupElement): Double = {
      val mixtureFrequency = mixture.get(Bases.basesToString(element.sequencedBases)).getOrElse(0.0)
      val probabilityCorrect = PhredUtils.phredToSuccessProbability(element.qualityScore) * element.read.alignmentLikelihood
      math.log10(mixtureFrequency * probabilityCorrect + (1-mixtureFrequency) * probabilityCorrect / 3.0 )
    }
    elements.map(logLikelihoodPileupElement _).sum
  }

  def findGermlineVariants(inputs: Seq[Input],
                           groupedInputs: GroupedInputs,
                           pileups: Seq[Pileup],
                          negativeLog10HeterozygousPrior: Double = 3,
                          negativeLog10HomozygousAlternatePrior: Double = 4,
                          negativeLog10CompoundAlternatePrior: Double = 8): Iterator[GermlineSmallVariant] = {

    val elements = groupedInputs.normalDNA.flatMap(pileups(_).elements)

    // We evaluate these models: (1) hom ref (2) het (3) hom alt (4) compound alt.
    // TODO: integrate strand bias into this likelihood.
    val ref = Bases.baseToString(pileups(0).referenceBase)
    val alleleReadCounts = elements.groupBy(element => Bases.basesToString(element.sequencedBases)).mapValues(_.size)
    val nonRefAlleles = alleleReadCounts.filterKeys(_ != ref).toSeq.sortBy(_._2).map(_._1)
    if (nonRefAlleles.isEmpty) {
      return Iterator.empty
    }

    val topAlt = nonRefAlleles(0)
    val secondAlt = if (nonRefAlleles.size > 1) nonRefAlleles(1) else "N"

    //val logMutationPrior = math.log10(mutationPrior)
    //val logNotMutationPrior = Seq(log10NotHeterozygousPrior, log10 math.log10(1 - mutationPrior)

    // TODO: get rid of this, and handle indels.
    if (topAlt.length != 1 || secondAlt.length != 1) {
      return Iterator.empty
    }

    // The following are not true posterior probabilities, but are proportional to posteriors.
    // Map from alleles to log probs.
    val logUnnormalizedPosteriors = Map(
      // Homozygous reference
      Seq(ref, ref) -> (logLikelihoodPileup(elements, Map(ref -> 1.0))),

      // Het
      Seq(topAlt) -> (logLikelihoodPileup(elements, Map(ref -> 0.5, topAlt -> 0.5)) - negativeLog10HeterozygousPrior),

      // Homozygous alt
      Seq(topAlt, topAlt) -> (logLikelihoodPileup(elements, Map(topAlt -> 1.0)) - negativeLog10HomozygousAlternatePrior),

      // Compound alt
      Seq(topAlt, secondAlt) ->
        (logLikelihoodPileup(elements, Map(topAlt -> 0.5, secondAlt -> 0.5)) - negativeLog10CompoundAlternatePrior)
    )
    val alts = logUnnormalizedPosteriors.toSeq.maxBy(_._2)._1
    if (alts == Seq(ref, ref))
      Iterator.empty  // no variants
    else
      Iterator(GermlineSmallVariant(pileups(0).referenceName, pileups(0).locus, ref, alts, logUnnormalizedPosteriors))
  }
}

