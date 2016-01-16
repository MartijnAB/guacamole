package org.hammerlab.guacamole.commands.jointcaller

import java.util

import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.variantcontext.{ Allele, GenotypeBuilder, VariantContext, VariantContextBuilder }
import htsjdk.variant.vcf._
import org.hammerlab.guacamole.DistributedUtil._
import org.hammerlab.guacamole.commands.jointcaller.SmallCharacterization.PileupStats
import org.hammerlab.guacamole.pileup.Pileup
import org.hammerlab.guacamole.reference.ReferenceBroadcast
import org.hammerlab.guacamole.{ Bases, ReadSet }

import scala.collection.{ mutable, JavaConversions }
import scala.collection.mutable.ArrayBuffer

case class SmallSomaticCharacterization(
  contig: String,
  position: Long,
  ref: String,
  topAlt: String,
  pooledCall: Boolean,
  individualCall: Boolean,
  germlineGenotypeQuality: Int,
  normalPercentError: Double,
  sampleMixtures: PerSample[Map[String, Double]],
  allelicDepthsPerSample: PerSample[Map[String, (Int, Int)]], // allele -> (total, positive strand)
  readNamesByAllele: PerSample[Map[String, Set[String]]])
    extends SmallCharacterization {

  def isVariant = pooledCall || individualCall

  def readNamesRef(): PerSample[Set[String]] = {
    readNamesByAllele.map(_.getOrElse(ref, Set.empty))
  }

  def readNamesAlt(): PerSample[Set[String]] = {
    readNamesByAllele.map(_.getOrElse(topAlt, Set.empty))
  }
}
object SmallSomaticCharacterization {

  case class Parameters(negativeLog10SomaticVariant: Double = 2,
                        minGermlineGenotypeQuality: Double = 10.0,
                        maxGermlineErrorRatePercent: Double = 4.0,
                        percentContaminationOfNormalWithTumor: Double = 0.0) {}
  val default = Parameters()

  def apply(germlineCharacterization: SmallGermlineCharacterization,
            normalDNAPileups: Seq[Pileup],
            tumorDNAPileups: Seq[Pileup],
            tumorRNAPileups: Seq[Pileup],
            parameters: Parameters = default): SmallSomaticCharacterization = {

    // We call a variant if we reject the null hypothesis in a likelihood ratio test for EITHER an individual
    // sample OR the pooled data.
    val individualDNASampleElements = tumorDNAPileups.map(_.elements.filter((_.read.alignmentQuality > 0)))
    val allDNAElements = individualDNASampleElements.flatten

    val ref = germlineCharacterization.ref

    val pooledStats = PileupStats(ref, allDNAElements)
    val sampleStats = individualDNASampleElements.map(elements => PileupStats(ref, elements))
    val topAlt = pooledStats.topAlt

    def unnormalizedPosteriors(stats: PileupStats): Map[Map[String, Double], Double] = {
      val topAltVaf = stats.alleleStats.get(topAlt)
        .map(_._3.size / stats.downsampledElements.size.toDouble).getOrElse(0.0)

      val singleAltMixture = Map(ref -> (1.0 - topAltVaf), topAlt -> topAltVaf)

      if (germlineCharacterization.isVariant) {
        Map(
          Map(ref -> 1.0) -> 0.0
        )
      } else {
        Map(
          Map(ref -> 1.0) -> stats.logLikelihoodPileup(germlineCharacterization.genotypeVafs),

          // Somatic
          singleAltMixture -> (stats.logLikelihoodPileup(singleAltMixture) - parameters.negativeLog10SomaticVariant))
      }
    }

    def call(posteriors: Map[Map[String, Double], Double]): Map[String, Double] = posteriors.maxBy(_._2)._1

    val normalStats = PileupStats(ref, normalDNAPileups.map(_.elements).flatten)
    val normalGenotypeDepth =
      Seq(germlineCharacterization.genotype._1, germlineCharacterization.genotype._2).distinct.map(
        allele => normalStats.allelicDepths.get(allele).map(_._1).getOrElse(0)).sum
    val normalTotalDepth = normalStats.allelicDepths.map(_._2._1).sum
    val normalErrors = normalTotalDepth -
      normalGenotypeDepth -
      (parameters.percentContaminationOfNormalWithTumor / 100.0
        * normalStats.allelicDepths.get(topAlt).map(_._1).getOrElse(0))
    val normalErrorRate = normalErrors.toDouble / math.max(normalTotalDepth, 1)

    val noCallGenotype = Map(ref -> 1.0)
    def noCallGenotypesPerSample = tumorDNAPileups.map(pileup => noCallGenotype)

    def makeResult(triggerPooled: Boolean, triggerIndividual: Boolean, perSampleGenotypes: PerSample[Map[String, Double]]): SmallSomaticCharacterization = {
      SmallSomaticCharacterization(
        tumorDNAPileups(0).referenceName,
        tumorDNAPileups(0).locus,
        ref,
        topAlt,
        triggerPooled,
        triggerIndividual,
        germlineCharacterization.genotypeQuality,
        normalErrorRate * 100.0,
        perSampleGenotypes,
        sampleStats.map(_.allelicDepths),
        sampleStats.map(_.readNamesByAllele)
      )
    }

    if (germlineCharacterization.isVariant) {
      // TODO: call LOH or reversion to normal.
      makeResult(false, false, noCallGenotypesPerSample)
    } else {
      val normalDNASupportsCalling = germlineCharacterization.genotypeQuality >= parameters.minGermlineGenotypeQuality &&
        normalErrorRate <= parameters.maxGermlineErrorRatePercent / 100.0

      if (!normalDNASupportsCalling) {
        // Not a high quality enough germline call to do a somatic call.
        makeResult(false, false, noCallGenotypesPerSample)
      } else {
        // Calling a regular somatic SNV.
        val pooledPosteriors = unnormalizedPosteriors(pooledStats)
        val perSamplePosteriors = sampleStats.map(unnormalizedPosteriors _)
        val perSampleGenotypes = perSamplePosteriors.map(call _)
        makeResult(call(pooledPosteriors) != noCallGenotype, perSampleGenotypes.exists(_ != noCallGenotype), perSampleGenotypes)
      }
    }
  }

  /**
   *
   * @param haplotype consecutive characterizations that form the same VCF-style variant
   * @param sampleNames
   * @param reference
   * @return
   */
  def makeHtsjdVariantContext(haplotype: Seq[SmallSomaticCharacterization],
                              sampleNames: PerSample[String],
                              reference: ReferenceBroadcast): VariantContext = {
    assume(haplotype.nonEmpty)
    val contig = haplotype.head.contig
    assume(haplotype.forall(item => item.contig == haplotype(0).contig))
    assume(haplotype.map(_.position) == haplotype.head.position.to(haplotype.last.position))

    val rawAltSequence = haplotype.map(_.topAlt).mkString
    val refPrepend = if (rawAltSequence.isEmpty) {
      // VCF requires non-empty alt alleles (e.g. a deletion C -> "" is coded as AC -> "A" where A is the preceding base)
      Bases.baseToString(reference.getContig(contig)(haplotype(0).position.toInt - 1))
    } else {
      ""
    }
    val altSequence = refPrepend + rawAltSequence
    val start_inclusive = haplotype.head.position - refPrepend.length
    val end_inclusive = haplotype.last.position
    assert(end_inclusive - start_inclusive + 1 == refPrepend.length + haplotype.length)

    val refSequence = Bases.basesToString(reference.getContig(contig)
      .slice(start_inclusive.toInt, end_inclusive.toInt + 1))
    assert(refPrepend + haplotype.map(_.ref).mkString == refSequence)

    val alleles = Seq(refSequence, altSequence).distinct

    def makeHtsjdkAllele(allele: String): Allele = Allele.create(allele, allele == refSequence)

    def range[K: Numeric](f: SmallSomaticCharacterization => K): String = {
      val values = haplotype.map(f(_))
      val max = values.max
      val min = values.min
      if (max == min) max.toString else "%s-%s".format(min.toString, max.toString)
    }

    val genotypes = sampleNames.zipWithIndex.map({
      case (name, sampleIndex) => {
        //val allelicDepths = items.map(_.allelicDepthsPerSample(sampleIndex).map(kv => refB(kv._1) -> kv._2))
        //val sampleMixture = range(_.sampleMixtures(sampleIndex))
        ///val sampleAlleles = sampleMixture.keys.toSeq.map(adjustAllele)

        // Are there any reads that imply this haplotype for this sample?
        val called = haplotype.map(_.readNamesAlt()(sampleIndex)).reduce(_.intersect(_)).nonEmpty

        val duplicatedSampleAlleles = if (called) Seq(refSequence, altSequence) else Seq(refSequence, refSequence)

        new GenotypeBuilder(name)
          .alleles(JavaConversions.seqAsJavaList(duplicatedSampleAlleles.map(makeHtsjdkAllele _)))
          //.AD(alleles.map(allelicDepths.getOrElse(_, (0, 0))).map(_._1).toArray)
          .attribute("GGQ", range(_.germlineGenotypeQuality))
          .attribute("NPE", range(_.normalPercentError))
          // .PL(genotypeLikelihoods.toArray)
          // .DP(depth)
          // .GQ(genotypeQuality.toInt)
          // .attribute("RL", logUnnormalizedPosteriors.map(
          //   pair => "%s/%s=%1.2f".format(pair._1._1, pair._1._2, pair._2)).mkString(" "))
          // .attribute("SBP",
          //   alleles.map(allele => "%s=%1.2f".format(allele, strandBias.getOrElse(allele, 0.0))).mkString(" "))
          /*
          .attribute("ADP",
            alleles.map(allele => {
              val totalAndPositive = allelicDepths.getOrElse(allele, (0, 0))
              "%s=%d/%d".format(allele, totalAndPositive._2, totalAndPositive._1)
            }).mkString(" "))
            */
          .make
      }
    })

    val numPooled = haplotype.count(_.pooledCall)
    val numIndividual = haplotype.count(_.individualCall)
    def triggerBracket(num: Int): String = {
      if (num < haplotype.length) "[%d/%d]".format(num, haplotype.length) else ""
    }
    val trigger = ((if (numPooled > 0) Seq("POOLED" + triggerBracket(numPooled)) else Seq.empty) ++
      (if (numIndividual > 0) Seq("INDIVIDUAL" + triggerBracket(numIndividual)) else Seq.empty)).mkString("+")

    new VariantContextBuilder()
      .chr(contig)
      .start(start_inclusive + 1) // one based based inclusive
      .stop(end_inclusive + 1)
      .genotypes(JavaConversions.seqAsJavaList(genotypes))
      .alleles(JavaConversions.seqAsJavaList(alleles.map(makeHtsjdkAllele _)))
      .attribute("TRIGGER", if (trigger.nonEmpty) trigger else "NONE")
      .make
  }

  def writeVcf(
    path: String,
    sampleNames: Seq[String],
    readSets: Seq[ReadSet],
    calls: Iterable[SmallSomaticCharacterization],
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
    headerLines.add(new VCFFormatHeaderLine("GGQ", 1, VCFHeaderLineType.Integer, "germline call genotype quality"))
    headerLines.add(new VCFFormatHeaderLine("NPE", 1, VCFHeaderLineType.Integer, "percent errors in normal"))

    // INFO
    headerLines.add(new VCFInfoHeaderLine("TRIGGER", 1, VCFHeaderLineType.String,
      "Which likelihood ratio test triggered call: POOLED and/or INDIVIDUAL"))

    val header = new VCFHeader(headerLines, JavaConversions.seqAsJavaList(sampleNames))
    header.setSequenceDictionary(readSets(0).sequenceDictionary.get.toSAMSequenceDictionary)
    header.setWriteCommandLine(true)
    header.setWriteEngineHeaders(true)
    header.addMetaDataLine(new VCFHeaderLine("Caller", "Guacamole"))
    writer.writeHeader(header)

    val groupedCalls = mutable.ArrayBuffer.empty[SmallSomaticCharacterization]

    def intersectSize(group1: PerSample[Set[String]], group2: PerSample[Set[String]]): Int = {
      group1.zip(group2).map(pair => pair._1.count(pair._2.contains _)).sum
    }

    def haplotypesWithEvidence(prefix: Seq[SmallSomaticCharacterization],
                               remaining: Seq[SmallSomaticCharacterization]): Seq[Seq[SmallSomaticCharacterization]] = {
      if (remaining.isEmpty) {
        Seq(prefix)
      } else {
        val minPosition = remaining.map(_.position).min
        val (potentialAdditions, rest) = remaining.partition(_.position == minPosition)
        potentialAdditions.flatMap(addition => {
          val (numSameHaplotype, numEndCurrentHaplotype, numStartNewHaplotype) = if (prefix.nonEmpty) {
            assert(minPosition == prefix.last.position + 1)
            (intersectSize(prefix.last.readNamesAlt, addition.readNamesAlt),
              intersectSize(prefix.last.readNamesAlt, addition.readNamesRef),
              intersectSize(prefix.last.readNamesRef, addition.readNamesAlt))
          } else {
            (0, 0, 1) // Start a new haplotype.
          }
          (if (numSameHaplotype > 0) haplotypesWithEvidence(prefix ++ Seq(addition), rest) else Seq.empty) ++
            (if (numEndCurrentHaplotype > 0) Seq(prefix) else Seq.empty) ++
            (if (numStartNewHaplotype > 0) haplotypesWithEvidence(Seq(addition), rest) else Seq.empty)
        })
      }
    }

    val variantContexts = ArrayBuffer.newBuilder[VariantContext]

    def processGroup(calls: Seq[SmallSomaticCharacterization]): Unit = {
      val haplotypes = haplotypesWithEvidence(Seq.empty, calls)
      variantContexts ++= haplotypes.map(makeHtsjdVariantContext(_, sampleNames, reference))
    }

    calls.foreach(call => {
      val groupWithPrevious = call.isVariant &&
        groupedCalls.lastOption.forall(
          prev => prev.contig == call.contig && call.position - prev.position <= 1)

      if (groupWithPrevious) {
        groupedCalls.append(call)
      } else {
        processGroup(groupedCalls)
        groupedCalls.clear()
        groupedCalls.append(call)
      }
    })
    if (groupedCalls.nonEmpty) {
      processGroup(groupedCalls)
      groupedCalls.clear()
    }

    variantContexts.result.sortBy(context => (context.getContig, context.getStart)).foreach(writer.add(_))
    writer.close()

  }
}