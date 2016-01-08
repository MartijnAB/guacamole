package org.hammerlab.guacamole.commands.jointcaller

import java.util

import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.variantcontext.{ GenotypeBuilder, VariantContextBuilder, Allele, VariantContext }
import htsjdk.variant.vcf._
import org.apache.spark.SparkContext
import org.bdgenomics.adam.util.PhredUtils
import org.hammerlab.guacamole.DistributedUtil.PerSample
import org.hammerlab.guacamole.commands.jointcaller.GenomicCharacterization.GermlineCharacterization.Parameters
import org.hammerlab.guacamole.commands.jointcaller.Inputs.{ GroupedInputs, Input }
import org.hammerlab.guacamole.commands.jointcaller.SomaticJoint.Arguments
import org.hammerlab.guacamole._
import org.hammerlab.guacamole.pileup.{ PileupElement, Pileup }
import org.hammerlab.guacamole.pileup.{ PileupElement, Pileup }
import org.hammerlab.guacamole.reads.Read
import org.hammerlab.guacamole.reference.ReferenceBroadcast
import sun.reflect.generics.reflectiveObjects.NotImplementedException

import scala.collection.mutable.ArrayBuffer
import scala.collection.{ mutable, JavaConversions }

trait GenomicCharacterization {
  def toHtsjdVariantContext(sampleNames: Seq[String]): VariantContext
}

object GenomicCharacterization {
  def logLikelihoodPileup(elements: Iterable[PileupElement], mixture: Map[String, Double]): Double = {
    def logLikelihoodPileupElement(element: PileupElement): Double = {
      if (element.read.alignmentQuality == 0) {
        0.0
      } else {
        val mixtureFrequency = mixture.get(Bases.basesToString(element.sequencedBases)).getOrElse(0.0)

        val probabilityCorrect =
          PhredUtils.phredToSuccessProbability(element.qualityScore) * element.read.alignmentLikelihood * .99

        val loglikelihood = math.log10(
          mixtureFrequency * probabilityCorrect + (1 - mixtureFrequency) * (1 - probabilityCorrect))
        loglikelihood
      }
    }
    val loglikelihoods = elements.map(logLikelihoodPileupElement _)
    loglikelihoods.sum
  }

  case class SomaticCharacterization(contig: String,
                                     position: Long,
                                     ref: String,
                                     topAlt: String,
                                     isVariant: Boolean,
                                     allelicDepthsPerSample: PerSample[Map[String, (Int, Int)]])  // allele -> (total, positive strand)
                                      extends GenomicCharacterization {

    def toHtsjdVariantContext(sampleNames: Seq[String]): VariantContext = {
      def makeHtsjdkAllele(allele: String): Allele = Allele.create(allele, allele == ref)

      val alleles = Seq(ref, topAlt).distinct

      val genotypes = sampleNames.zip(allelicDepthsPerSample).map({case (name, allelicDepths) => {
        new GenotypeBuilder(name)
          .alleles(JavaConversions.seqAsJavaList(Seq(ref, topAlt).map(makeHtsjdkAllele _)))
          .AD(alleles.map(allelicDepths.getOrElse(_, (0, 0))).map(_._1).toArray)
            // .PL(genotypeLikelihoods.toArray)
            // .DP(depth)
            // .GQ(genotypeQuality.toInt)
            // .attribute("RL", logUnnormalizedPosteriors.map(
            //   pair => "%s/%s=%1.2f".format(pair._1._1, pair._1._2, pair._2)).mkString(" "))
            // .attribute("SBP",
            //   alleles.map(allele => "%s=%1.2f".format(allele, strandBias.getOrElse(allele, 0.0))).mkString(" "))
          .attribute("ADP",
            alleles.map(allele => {
              val totalAndPositive = allelicDepths.getOrElse(allele, (0, 0))
              "%s=%d/%d".format(allele, totalAndPositive._2, totalAndPositive._1)
            }).mkString(" "))
            .make
      }})

      new VariantContextBuilder()
        .chr(contig)
        .start(position + 1) // one based based inclusive
        .stop(position + 1 + math.max(ref.length - 1, 0))
        .genotypes(JavaConversions.seqAsJavaList(genotypes))
        .alleles(JavaConversions.seqAsJavaList(alleles.map(makeHtsjdkAllele _)))
        .make
    }

  }
  object SomaticCharacterization {

    case class Parameters(negativeLog10SomaticVariant: Double = 2,
                           minGermlineGenotypeQuality: Double = 10.0)
    {}
    val default = Parameters()

    def apply(germlineCharacterization: GermlineCharacterization,
              tumorDNAPileups: Seq[Pileup],
              tumorRNAPileups: Seq[Pileup],
              parameters: Parameters = default): SomaticCharacterization = {

      // We call a variant if we reject the null hypothesis in a likelhood ratio test for EITHER an individual
      // sample OR the pooled data.

      val individualDNASampleElements = tumorDNAPileups.map(_.elements.filter((_.read.alignmentQuality > 0)))
      val allDNAElements = individualDNASampleElements.flatten

      val ref = germlineCharacterization.ref

      var topAltPooled: Option[String] = None

      val (pooledDepths, pooledCall) :: sampleDepthsAndCalls = (Seq(allDNAElements) ++ individualDNASampleElements).map(elements => {

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

        // First iteratorion of loop we're processing the pooled data, so we store topAlt here.
        val topAlt = topAltPooled.getOrElse(nonRefAlleles.headOption.getOrElse("N"))
        topAltPooled = Some(topAlt)

        val topAltVaf = allelicDepthsAndSubsampled.get(
          topAlt).map(_._3.size / downsampledElements.size.toDouble).getOrElse(0.0)

        val singleAltMixture = Map(ref -> (1.0 - topAltVaf), topAlt -> topAltVaf)

        val logUnnormalizedPosteriors = if (germlineCharacterization.isVariant) {
          Map(
            Map(ref -> 1.0) -> 0.0
          )
        } else {
          Map(
            Map(ref -> 1.0) -> logLikelihoodPileup(downsampledElements, germlineCharacterization.genotypeVafs),

            // Somatic
            singleAltMixture -> (logLikelihoodPileup(downsampledElements, singleAltMixture)
              - parameters.negativeLog10SomaticVariant))
        }
        (allelicDepths, logUnnormalizedPosteriors.toSeq.maxBy(_._2)._1)
      })

      val noCallGenotype = Map(ref -> 1.0)

      SomaticCharacterization(
        tumorDNAPileups(0).referenceName,
        tumorDNAPileups(0).locus,
        ref,
        topAltPooled.get,
        pooledCall != noCallGenotype || sampleDepthsAndCalls.exists(_._2 != noCallGenotype),
        sampleDepthsAndCalls.map(_._1))
    }

    def writeVcf(
                  path: String,
                  sampleNames: Seq[String],
                  readSets: Seq[ReadSet],
                  calls: Iterable[SomaticCharacterization],
                  reference: ReferenceBroadcast) = {
      // TODO

      throw new NotImplementedException()
    }
  }

  case class GermlineCharacterization(
    contig: String,
    position: Long,
    ref: String,
    topAlt: String,
    secondAlt: String,
    genotype: (String, String),
    allelicDepths: Map[String, (Int, Int)], // allele -> (total, positive strand)
    depth: Int,
    logUnnormalizedPosteriors: Map[(String, String), Double])
      extends GenomicCharacterization {

    def topAltAllelicFraction =
      Seq(genotype._1, genotype._2).count(_ == topAlt) * 0.5

    def genotypeVafs =
      if (genotype._1 == genotype._2) Map(genotype._1 -> 1.0) else Map(genotype._1 -> 0.5, genotype._2 -> 0.5)

    def isVariant = genotype != (ref, ref) && topAlt != "N"

    def possibleGenotypes = Seq((ref, ref), (ref, topAlt), (topAlt, topAlt)) ++
      (if (secondAlt != "N") Seq((topAlt, secondAlt)) else Seq.empty)

    def genotypeLikelihoods = GermlineCharacterization.logPosteriorsToNormalizedPhred(
      possibleGenotypes.map(
        value => logUnnormalizedPosteriors.getOrElse(value, Double.NegativeInfinity))).map(_.toInt)

    def genotypeQuality = genotypeLikelihoods.sorted.seq(1)

    def toHtsjdVariantContext(sampleNames: Seq[String]): VariantContext = {
      def makeHtsjdkAllele(allele: String): Allele = Allele.create(allele, allele == ref)

      assume(sampleNames.size == 1)
      val sampleName = sampleNames(0)

      val alleles = Seq(ref, genotype._1, genotype._2).distinct

      // If there is a true second alt with evidence, the PL field will have 4 values (hom ref, het, hom alt, compound
      // alt), otherwise just 3. GATK always gives the 3-genotype version.

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
                          negativeLog10CompoundAlternatePrior: Double = 4) {}
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
      sampleNames: Seq[String],
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

      val header = new VCFHeader(headerLines, JavaConversions.seqAsJavaList(sampleNames))
      header.setSequenceDictionary(readSets(0).sequenceDictionary.get.toSAMSequenceDictionary)
      header.setWriteCommandLine(true)
      header.setWriteEngineHeaders(true)
      header.addMetaDataLine(new VCFHeaderLine("Caller", "Guacamole"))
      writer.writeHeader(header)

      val normalizedCalls = normalizeDeletions(reference, calls)

      var prevChr = ""
      var prevStart = -1
      normalizedCalls.foreach(call => {
        val variantContext = call.toHtsjdVariantContext(sampleNames)
        if (prevChr == variantContext.getContig) {
          assert(variantContext.getStart >= prevStart,
            "Out of order: expected %d >= %d".format(variantContext.getStart, prevStart))
        }
        prevChr = variantContext.getContig
        prevStart = variantContext.getStart
        writer.add(variantContext)
      })
      writer.close()
    }
  }
}
