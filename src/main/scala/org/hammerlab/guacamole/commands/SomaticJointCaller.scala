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
import org.apache.commons.math3.util.ArithmeticUtils
import org.apache.http.client.utils.URLEncodedUtils
import org.apache.spark.SparkContext
import org.bdgenomics.adam.util.PhredUtils
import org.hammerlab.guacamole.Common.Arguments.NoSequenceDictionary
import org.hammerlab.guacamole.pileup._
import org.hammerlab.guacamole.reads._
import org.hammerlab.guacamole._
import org.hammerlab.guacamole.reference.ReferenceBroadcast
import org.kohsuke.args4j.{ Option => Args4jOption, Argument }

import scala.collection.{mutable, JavaConversions}
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

    @Args4jOption(name = "--reference-fasta", required = false, usage = "Local path to a reference FASTA file")
    var referenceFastaPath: String = null

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
      if (element.read.alignmentQuality == 0) {
        0.0
      } else {
        val mixtureFrequency = mixture.get(Bases.basesToString(element.sequencedBases)).getOrElse(0.0)

        val probabilityCorrect =
          PhredUtils.phredToSuccessProbability(element.qualityScore) * element.read.alignmentLikelihood

        val loglikelihood = math.log10(
          mixtureFrequency * probabilityCorrect + (1 - mixtureFrequency) * (1 - probabilityCorrect))
        loglikelihood
      }
    }
    val loglikelihoods = elements.map(logLikelihoodPileupElement _)
    loglikelihoods.sum
  }

  case class GermlineCharacterization(
                                     contig: String,
                                     position: Long,
                                     ref: String,
                                     topAlt: String,
                                     secondAlt: String,
                                     genotype: (String, String),
                                     allelicDepths: Map[String, (Int, Int)],  // allele -> (total, positive strand)
                                     depth: Int,
                                     logUnnormalizedPosteriors: Map[(String, String), Double]) {

    def topAltAllelicFraction =
      Seq(genotype._1, genotype._2).count(_ == topAlt) * 0.5

    def isVariant = genotype != (ref, ref) && topAlt != "N"

    def toHtsjdVariantContext(sampleName: String): VariantContext = {
      def makeHtsjdkAllele(allele: String): Allele = Allele.create(allele, allele == ref)

      val alleles = Seq(ref, genotype._1, genotype._2).distinct

      // If there is a true second alt with evidence, the PL field will have 4 values (hom ref, het, hom alt, compound
      // alt), otherwise just 3. GATK always gives the 3-genotype version.
      val possibleGenotypes = Seq((ref, ref), (ref, topAlt), (topAlt, topAlt)) ++
        (if (secondAlt != "N") Seq((topAlt, secondAlt)) else Seq.empty)

      val genotypeLikelihoods = GermlineCharacterization.logPosteriorsToNormalizedPhred(
        possibleGenotypes.map(
          value => logUnnormalizedPosteriors.getOrElse(value, Double.NegativeInfinity))).map(_.toInt)

      val genotypeQuality = genotypeLikelihoods.sorted.seq(1)

      assert(genotype._1.nonEmpty)
      assert(genotype._2.nonEmpty)

      new VariantContextBuilder()
        .chr(contig)
        .start(position + 1) // one based based inclusive
        .stop(position + 1 + math.max(ref.length - 1, 0))
        .genotypes(
          new GenotypeBuilder(sampleName)
            .alleles(JavaConversions.seqAsJavaList(Seq(genotype._1, genotype._2).map(makeHtsjdkAllele _)))
            .AD(alleles.map(allelicDepths.getOrElse(_, (0, 0))).map(_._1).toArray)
            .PL(genotypeLikelihoods.toArray)
            .DP(depth)
            .GQ(genotypeQuality.toInt)
            .attribute("RL", logUnnormalizedPosteriors.map(
              pair => "%s/%s=%1.2f".format(pair._1._1, pair._1._2, pair._2)).mkString(" "))
           // .attribute("SBP",
           //   alleles.map(allele => "%s=%1.2f".format(allele, strandBias.getOrElse(allele, 0.0))).mkString(" "))
            .attribute("ADP",
              alleles.map(allele => {
                val totalAndPositive = allelicDepths.getOrElse(allele, (0, 0))
                "%s=%d/%d".format(allele, totalAndPositive._2, totalAndPositive._1)
              }).mkString(" "))
            .make)
        .alleles(JavaConversions.seqAsJavaList(alleles.map(makeHtsjdkAllele _)))
        .make
    }
  }

  object GermlineCharacterization {
    case class Parameters(negativeLog10HeterozygousPrior: Double = 2,
                          negativeLog10HomozygousAlternatePrior: Double = 2,
                          negativeLog10CompoundAlternatePrior: Double = 4)
    {}
    val default = Parameters()

    def apply(normalDNAPileups: Seq[Pileup], parameters: Parameters = default): GermlineCharacterization = {

      val elements = normalDNAPileups.flatMap(_.elements).filter(_.read.alignmentQuality > 0)
      val depth = elements.length
      val ref = Bases.baseToString(normalDNAPileups(0).referenceBase).toUpperCase

      // We evaluate these models: (1) hom ref (2) het (3) hom alt (4) compound alt.
      // Map from allele -> (num total supporting reads, num positive strand supporting reads, downsampled elements)
      val allelicDepthsAndSubsampled = elements.sortBy(_.qualityScore * -1).groupBy(element => {
        Bases.basesToString(element.sequencedBases).toUpperCase
      }).mapValues(items => {
        val positive = items.filter(_.read.isPositiveStrand)
        val negative = items.filter(!_.read.isPositiveStrand)
        val numToTake = Math.max(Math.min(positive.size, negative.size), 3)
        val adjustedItems = positive.take(numToTake) ++ negative.take(numToTake)
        (items.size, positive.size, adjustedItems)
      })
      val allelicDepths = allelicDepthsAndSubsampled.mapValues(item => (item._1, item._2))
      val downsampledElements = allelicDepthsAndSubsampled.mapValues(_._3).values.flatten

      val nonRefAlleles = allelicDepths.filterKeys(_ != ref).toSeq.sortBy(_._2._1 * -1).map(_._1)

      val topAlt = nonRefAlleles.headOption.getOrElse("N")
      val secondAlt = if (nonRefAlleles.size > 1) nonRefAlleles(1) else "N"

      // The following are not true posterior probabilities, but are proportional to posteriors.
      // Map from alleles to log probs.
      val logUnnormalizedPosteriors = Map(
        // Homozygous reference.
        (ref, ref) -> (
          logLikelihoodPileup(downsampledElements, Map(ref -> 1.0))),

        // Het
        (ref, topAlt) -> (
          logLikelihoodPileup(downsampledElements, Map(ref -> 0.5, topAlt -> 0.5))
            - parameters.negativeLog10HeterozygousPrior),

        // Homozygous alt
        (topAlt, topAlt) -> (
          logLikelihoodPileup(downsampledElements, Map(topAlt -> 1.0))
            - parameters.negativeLog10HomozygousAlternatePrior),

        // Compound alt
        (topAlt, secondAlt) -> (
          logLikelihoodPileup(downsampledElements, Map(topAlt -> 0.5, secondAlt -> 0.5))
            - parameters.negativeLog10CompoundAlternatePrior)
      )
      val genotype = logUnnormalizedPosteriors.toSeq.maxBy(_._2)._1
      //if (genotype != (ref, ref)) {
      //  println("Adjusted %d -> %d".format(elements.size, downsampledElements.size))
      //}
      GermlineCharacterization(
        normalDNAPileups(0).referenceName,
        normalDNAPileups(0).locus,
        ref,
        topAlt,
        secondAlt,
        genotype,
        allelicDepths,
        depth,
        logUnnormalizedPosteriors)
    }

    def logPosteriorsToNormalizedPhred(log10Probabilities: Seq[Double], log10Precision: Double = -16): Seq[Double] = {
      val max = log10Probabilities.max
      val rescaled = log10Probabilities.map(_ - max)
      rescaled.map(_ * -10)
    }

    def normalizeDeletions(reference: ReferenceBroadcast, calls: Iterable[GermlineCharacterization]): Iterator[GermlineCharacterization] = {
      def makeDeletion(deletionCalls: Seq[GermlineCharacterization]): GermlineCharacterization = {
        val start = deletionCalls.head
        val end = deletionCalls.last
        assume(start.contig == end.contig)
        assume(end.position >= start.position)
        assume(start.topAlt.isEmpty)
        assume(end.topAlt.isEmpty)
        assume(deletionCalls.forall(_.topAltAllelicFraction == start.topAltAllelicFraction))

        val het = start.topAltAllelicFraction == 0.5

        val refSequence = Bases.basesToString(
          reference.getReferenceSequence(start.contig, start.position.toInt - 1, end.position.toInt + 1)).toUpperCase

        def convertAllele(allele: String) =
          refSequence(0) + allele

        GermlineCharacterization(
          start.contig,
          start.position - 1,
          refSequence,
          convertAllele(start.topAlt),
          convertAllele(start.secondAlt),
          if (het) (refSequence, convertAllele("")) else (convertAllele(""), convertAllele("")),
          start.allelicDepths.map(pair => convertAllele(pair._1) -> pair._2),
          start.depth,
          start.logUnnormalizedPosteriors.map(
            pair => (convertAllele(pair._1._1), convertAllele(pair._1._2)) -> pair._2))

      }

      val results = mutable.ArrayStack[GermlineCharacterization]()
      val currentDeletion = ArrayBuffer.empty[GermlineCharacterization]

      def maybeEmitDeletion(): Unit = {
        if (currentDeletion.nonEmpty) {
          // Ending a deletion
          val fullDeletion = makeDeletion(currentDeletion)
          if (results.headOption.exists(_.position == fullDeletion.position)) {
            // Combining with an existing call
            // TODO: someday, we may want to attempt to merge SNV and indel calls here.
            // For now, we just discard the non-indel call.
            results.pop()
          }
          results += fullDeletion
          currentDeletion.clear()
        }
      }

      calls.foreach(call => {
        val deletionContinuation = call.topAlt.isEmpty &&
            currentDeletion.lastOption.exists(last =>
              last.position == call.position - 1 &&
              last.contig == call.contig &&
              last.topAltAllelicFraction == call.topAltAllelicFraction)

        if (deletionContinuation) {
          // Continuing a deletion.
          currentDeletion += call
        } else {
          maybeEmitDeletion()
        }

        if (call.genotype == (call.topAlt, call.secondAlt) && call.secondAlt.isEmpty) {
          // Compound variant where the second alternate is a deletion. For now we skip these.
        } else if (call.topAlt.nonEmpty) {
          // A simple non-deletion call.
          results += call
        } else if (!deletionContinuation) {
          // Starting a deletion
          assert(currentDeletion.isEmpty)
          currentDeletion += call
        }
      })
      maybeEmitDeletion()
      results.foreach(result => {
        assert(result.genotype._1.nonEmpty)
        assert(result.genotype._2.nonEmpty)
      })
      results.reverseIterator
    }

    def writeVcf(
                  path: String,
                  sampleName: String,
                  readSets: Seq[ReadSet],
                  calls: Iterable[GermlineCharacterization],
                  reference: ReferenceBroadcast) = {
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
      headerLines.add(new VCFFormatHeaderLine("GQ", 1, VCFHeaderLineType.Integer, "Genotype quality"))

      // nonstandard
      headerLines.add(new VCFFormatHeaderLine("RL", 1, VCFHeaderLineType.String, "Unnormalized log10 genotype posteriors"))
      headerLines.add(new VCFFormatHeaderLine("SBP", 1, VCFHeaderLineType.String, "log10 strand bias p-values"))
      headerLines.add(new VCFFormatHeaderLine("ADP", 1, VCFHeaderLineType.String, "allelic depth as num postiive strand / num total"))

      val header = new VCFHeader(headerLines, JavaConversions.seqAsJavaList(Seq(sampleName)))
      header.setSequenceDictionary(readSets(0).sequenceDictionary.get.toSAMSequenceDictionary)
      header.setWriteCommandLine(true)
      header.setWriteEngineHeaders(true)
      header.addMetaDataLine(new VCFHeaderLine("Caller", "Guacamole"))
      writer.writeHeader(header)

      val normalizedCalls = normalizeDeletions(reference, calls)

      var prevChr = ""
      var prevStart = -1
      normalizedCalls.foreach(call => {
        val variantContext = call.toHtsjdVariantContext(sampleName)
        if (prevChr == variantContext.getChr) {
          assert(variantContext.getStart >= prevStart,
            "Out of order: expected %d >= %d".format(variantContext.getStart, prevStart))
        }
        prevChr = variantContext.getChr
        prevStart = variantContext.getStart
        writer.add(variantContext)
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

      val reference = Option(args.referenceFastaPath).map(ReferenceBroadcast(_, sc))
      if (reference.isEmpty) {
        throw new IllegalArgumentException("Reference fasta required")
      }

      val readSets = inputs.zipWithIndex.map({
        case (input, index) => ReadSet(
          sc,
          input.path,
          false,
          Read.InputFilters.empty,
          token = index,
          contigLengthsFromDictionary = !args.noSequenceDictionary,
          referenceGenome = reference)
      })

      assert(readSets.forall(_.sequenceDictionary == readSets(0).sequenceDictionary),
        "Samples have different sequence dictionaries: %s."
          .format(readSets.map(_.sequenceDictionary.toString).mkString("\n")))

      val loci = Common.lociFromArguments(args, readSets(0))
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

      val broadcastForceCallLoci = sc.broadcast(forceCallLoci)
      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(args, loci, readSets.map(_.mappedReads): _*)
      val groupedInputs = GroupedInputs(inputs)

      // TODO: we currently are re-shuffling on every pileupFlatMap call.
      // Call germline small variants.
      val germlineCalls = DistributedUtil.pileupFlatMapMultipleRDDs(
        groupedInputs.normalDNA.map(readSets(_).mappedReads),
        lociPartitions,
        false, // skip empty
        pileups => {
          val characterization = GermlineCharacterization(groupedInputs.normalDNA.map(pileups(_)))
          val emit = (characterization.isVariant
            || broadcastForceCallLoci.value.onContig(pileups(0).referenceName).contains(pileups(0).locus))
          if (emit) Iterator(characterization) else Iterator.empty
        }, referenceGenome = reference).collect

      Common.progress("Called %,d germline variants.".format(germlineCalls.length))

      if (args.outSmallGermlineVariants.nonEmpty) {
        Common.progress("Writing germline variants.")
        GermlineCharacterization.writeVcf(
          args.outSmallGermlineVariants, args.normalSampleName, readSets, germlineCalls, reference.get)
      }
      Common.progress("Wrote %,d calls to %s".format(germlineCalls.length, args.outSmallGermlineVariants))
    }
  }

  def combinedPileup(pileups: Seq[Pileup]) = {
    val elements = pileups.flatMap(_.elements)
    Pileup(pileups(0).referenceName, pileups(0).locus, pileups(0).referenceBase, elements)
  }
}
