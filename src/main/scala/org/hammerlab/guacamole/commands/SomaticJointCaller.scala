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

import java.util

import htsjdk.samtools.util.Locus
import htsjdk.variant.variantcontext.{ GenotypeBuilder, Allele, VariantContextBuilder, VariantContext }
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf._
import org.apache.http.client.utils.URLEncodedUtils
import org.apache.spark.SparkContext
import org.bdgenomics.adam.util.PhredUtils
import org.hammerlab.guacamole.Common.Arguments.NoSequenceDictionary
import org.hammerlab.guacamole.pileup._
import org.hammerlab.guacamole.reads._
import org.hammerlab.guacamole._
import org.kohsuke.args4j.{ Option => Args4jOption, Argument }

import scala.collection.JavaConversions
import scala.collection.mutable.ArrayBuffer

object SomaticJoint {
  class Arguments extends DistributedUtil.Arguments with NoSequenceDictionary {

    @Argument(required = true, multiValued = true,
      usage = "FILE1 FILE2 FILE3")
    var inputs: Array[String] = Array.empty

    @Args4jOption(name = "--force-call-loci-from-file", usage = "Always call the given sites")
    var forceCallLociFromFile: String = ""

    @Args4jOption(name = "--force-call-loci", usage = "Always call the given sites")
    var forceCallLoci: String = ""

    @Args4jOption(name = "--normal-sample-name", usage = "Sample name to use for VCF output of germline calls")
    var normalSampleName: String = "NORMAL"

    @Args4jOption(name = "--out-small-germline-variants", usage = "Output path. Default: no output")
    var outSmallGermlineVariants: String = ""

    @Args4jOption(name = "--out-small-somatic-variants", usage = "Output path. Default: no outupt")
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

  def logLikelihoodPileup(elements: Iterable[PileupElement], mixture: Map[String, Double]): Double = {
    def logLikelihoodPileupElement(element: PileupElement): Double = {
      val mixtureFrequency = mixture.get(Bases.basesToString(element.sequencedBases)).getOrElse(0.0)
      val probabilityCorrect =
        PhredUtils.phredToSuccessProbability(element.qualityScore) * element.read.alignmentLikelihood
      math.log10(mixtureFrequency * probabilityCorrect + (1 - mixtureFrequency) * probabilityCorrect / 3.0)
    }
    elements.map(logLikelihoodPileupElement _).sum
  }

  class GermlineSmallVariantCall(negativeLog10HeterozygousPrior: Double = 2,
                                 negativeLog10HomozygousAlternatePrior: Double = 4,
                                 negativeLog10CompoundAlternatePrior: Double = 8) {

    // Essential information.
    var contig = ""
    var position = -1L
    var ref = ""
    var topAlt = ""
    var secondAlt = ""
    var genotype = ("N", "N")

    // Additional diagnostic info.
    var allelicDepths = Map[String, Int]()
    var depth = -1
    var logUnnormalizedPosteriors = Map[(String, String), Double]()

    def call(normalDNAPileups: Seq[Pileup]): Boolean = {
      val elements = normalDNAPileups.flatMap(_.elements)
      depth = elements.length

      contig = normalDNAPileups(0).referenceName
      position = normalDNAPileups(0).locus
      ref = Bases.baseToString(normalDNAPileups(0).referenceBase)

      // We evaluate these models: (1) hom ref (2) het (3) hom alt (4) compound alt.
      // TODO: integrate strand bias into this likelihood.
      allelicDepths = elements.groupBy(element => Bases.basesToString(element.sequencedBases)).mapValues(_.size)
      val nonRefAlleles = allelicDepths.filterKeys(_ != ref).toSeq.sortBy(_._2 * -1).map(_._1)

      topAlt = nonRefAlleles.headOption.getOrElse("N")
      secondAlt = if (nonRefAlleles.size > 1) nonRefAlleles(1) else "N"

      //val logMutationPrior = math.log10(mutationPrior)
      //val logNotMutationPrior = Seq(log10NotHeterozygousPrior, log10 math.log10(1 - mutationPrior)

      // TODO: get rid of this, and handle indels.
      if (topAlt.length != 1 || secondAlt.length != 1) {
        return false
      }

      // The following are not true posterior probabilities, but are proportional to posteriors.
      // Map from alleles to log probs.
      logUnnormalizedPosteriors = Map(
        // Homozygous reference
        (ref, ref) -> (logLikelihoodPileup(elements, Map(ref -> 1.0))),

        // Het
        (ref, topAlt) -> (logLikelihoodPileup(elements, Map(ref -> 0.5, topAlt -> 0.5)) - negativeLog10HeterozygousPrior),

        // Homozygous alt
        (topAlt, topAlt) -> (logLikelihoodPileup(elements, Map(topAlt -> 1.0)) - negativeLog10HomozygousAlternatePrior),

        // Compound alt
        (topAlt, secondAlt) ->
          (logLikelihoodPileup(elements, Map(topAlt -> 0.5, secondAlt -> 0.5)) - negativeLog10CompoundAlternatePrior)
      )
      genotype = logUnnormalizedPosteriors.toSeq.maxBy(_._2)._1
      genotype != (ref, ref)
    }

    def toHtsjdVariantContext(sampleName: String): VariantContext = {
      def makeHtsjdkAllele(allele: String): Allele = Allele.create(allele, allele == ref)

      val alleles = Seq(ref, genotype._1, genotype._2).distinct

      // If there is a true second alt with evidence, the PL field will have 4 values (hom ref, het, hom alt, compound
      // alt), otherwise just 3. The 4-allele version seems to be different than what e.g. the GATK does, but is
      // more useful.
      val possibleGenotypes = Seq((ref, ref), (ref, topAlt), (topAlt, topAlt)) ++
        (if (secondAlt != "N") Seq((topAlt, secondAlt)) else Seq.empty)

      new VariantContextBuilder()
        .chr(contig)
        .start(position + 1) // one based based inclusive
        .stop(position + 1 + math.max(ref.length - 1, 0))
        .genotypes(
          new GenotypeBuilder(sampleName)
            .alleles(JavaConversions.seqAsJavaList(Seq(genotype._1, genotype._2).map(makeHtsjdkAllele _)))
            .AD(alleles.map(allelicDepths.getOrElse(_, 0)).toArray)
            .PL(GermlineSmallVariantCall.unnormalizedLogProbsToNormalizedPhred(
              possibleGenotypes.map(
                value => logUnnormalizedPosteriors.getOrElse(value, Double.NegativeInfinity))).toArray)
            .DP(depth)
            .attribute("RL", logUnnormalizedPosteriors.map(
              pair => "%s/%s=%1.2f".format(pair._1._1, pair._1._2, pair._2)).mkString(" "))
            .make)
        .alleles(JavaConversions.seqAsJavaList(alleles.map(makeHtsjdkAllele _)))
        .make
    }
  }

  object GermlineSmallVariantCall {
    def unnormalizedLogProbsToNormalizedPhred(log10Probabilities: Seq[Double], log10Precision: Double = -16): Seq[Double] = {
      // See: http://stats.stackexchange.com/questions/66616/converting-normalizing-very-small-likelihood-values-to-probability
      // val threshold = log10Precision - math.log10(log10Probabilities.length.toDouble)
      val max = log10Probabilities.max
      val rescaled = log10Probabilities.map(_ - max)
      val exponentials = rescaled.map(math.pow(10, _))
      val sum = exponentials.sum
      exponentials.map(value => math.log10(1 - (value / sum)) * -10.0)
    }

    def writeVcf(path: String, sampleName: String, readSets: Seq[ReadSet], calls: Iterable[GermlineSmallVariantCall]) = {
      val writer = new VariantContextWriterBuilder()
        .setOutputFile(path)
        .setReferenceDictionary(readSets(0).sequenceDictionary.get.toSAMSequenceDictionary)
        .build
      val headerLines = new util.HashSet[VCFHeaderLine]()
      headerLines.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"))
      headerLines.add(new VCFFormatHeaderLine("AD", VCFHeaderLineCount.R, VCFHeaderLineType.Integer,
        "Allelic depths for the ref and alt alleles"))
      headerLines.add(new VCFFormatHeaderLine("PL", VCFHeaderLineCount.G, VCFHeaderLineType.Integer,
        "Phred scaled genotype likelihoods"))
      headerLines.add(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Total depth"))

      // nonstandard
      headerLines.add(new VCFFormatHeaderLine("RL", 1, VCFHeaderLineType.String, "Unnormalized log10 genotype posteriors"))

      val header = new VCFHeader(headerLines, JavaConversions.seqAsJavaList(Seq(sampleName)))
      header.setSequenceDictionary(readSets(0).sequenceDictionary.get.toSAMSequenceDictionary)
      header.setWriteCommandLine(true)
      header.setWriteEngineHeaders(true)
      header.addMetaDataLine(new VCFHeaderLine("Caller", "Guacamole"))
      writer.writeHeader(header)

      calls.foreach(call => {
        writer.add(call.toHtsjdVariantContext(sampleName))
      })
      writer.close()
    }
  }

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

      val loci = Common.lociFromArguments(args)
      val forceCallLoci = if (args.forceCallLoci.nonEmpty || args.forceCallLociFromFile.nonEmpty) {
        Common.loci(args.forceCallLoci, args.forceCallLociFromFile, readSets(0))
      } else {
        LociSet.empty
      }

      if (forceCallLoci.nonEmpty) {
        Common.progress("Force calling %,d loci across %,d contig(s): %s".format(
          forceCallLoci.count,
          forceCallLoci.contigs.length,
          forceCallLoci.truncatedString()))
      }

      readSets(0).source

      val broadcastForceCallLoci = sc.broadcast(forceCallLoci)
      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(
        args, loci.result(readSets(0).contigLengths), readSets.map(_.mappedReads): _*)
      val groupedInputs = GroupedInputs(inputs)

      // TODO: we currently are re-shuffling on every pileupFlatMap call.
      // Call germline small variants.
      val germlineCalls = DistributedUtil.pileupFlatMapMultipleRDDs(
        groupedInputs.normalDNA.map(readSets(_).mappedReads),
        lociPartitions,
        false, // skip empty
        pileups => {
          val maybeCall = new GermlineSmallVariantCall()
          val called = maybeCall.call(groupedInputs.normalDNA.map(pileups(_)))
          if (called || broadcastForceCallLoci.value.onContig(pileups(0).referenceName).contains(pileups(0).locus)) {
            Iterator(maybeCall)
          } else {
            Iterator.empty
          }
        }).collect

      Common.progress("Called %,d germline variants.".format(germlineCalls.length))

      if (args.outSmallGermlineVariants.nonEmpty) {
        Common.progress("Writing germline variants.")
        GermlineSmallVariantCall.writeVcf(args.outSmallGermlineVariants, args.normalSampleName, readSets, germlineCalls)
      }
      Common.progress("Wrote %,d calls to %s".format(germlineCalls.length, args.outSmallGermlineVariants))
    }
  }

  def combinedPileup(pileups: Seq[Pileup]) = {
    val elements = pileups.flatMap(_.elements)
    Pileup(pileups(0).referenceName, pileups(0).locus, pileups(0).referenceBase, elements)
  }
}
