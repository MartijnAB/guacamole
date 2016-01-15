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

import scala.collection.JavaConversions

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
  allelicDepthsPerSample: PerSample[Map[String, (Int, Int)]],
  readSupportRef: PerSample[Set[String]],
  readSupportAlt: PerSample[Set[String]]) // allele -> (total, positive strand)
    extends SmallCharacterization {

  def isVariant = pooledCall || individualCall

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
        sampleStats.map(_.readsSupportingAllele(ref)),
        sampleStats.map(_.readsSupportingAllele(topAlt))
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

  def makeHtsjdVariantContext(item: SmallSomaticCharacterization, sampleNames: Seq[String], reference: ReferenceBroadcast): VariantContext = {
    val (adjustedPosition, adjustAllele) = if (item.topAlt.nonEmpty) {
      (item.position, (allele: String) => allele)
    } else {
      val refSeqeunce = Bases.basesToString(
        reference.getReferenceSequence(item.contig, item.position.toInt - 1, item.position.toInt + 1)).toUpperCase
      (item.position - 1, (allele: String) => refSeqeunce(0) + allele)
    }

    val alleles = Seq(adjustAllele(item.ref), adjustAllele(item.topAlt)).distinct

    def makeHtsjdkAllele(allele: String): Allele = Allele.create(allele, allele == adjustAllele(item.ref))

    val genotypes = sampleNames.zipWithIndex.map({
      case (name, sampleIndex) => {
        val allelicDepths = item.allelicDepthsPerSample(sampleIndex).map(kv => (adjustAllele(kv._1) -> kv._2))
        val sampleMixture = item.sampleMixtures(sampleIndex)
        val sampleAlleles = sampleMixture.keys.toSeq.map(adjustAllele)
        val duplicatedSampleAlleles = if (sampleAlleles.length == 1)
          Seq(sampleAlleles(0), sampleAlleles(0))
        else
          sampleAlleles

        new GenotypeBuilder(name)
          .alleles(JavaConversions.seqAsJavaList(duplicatedSampleAlleles.map(makeHtsjdkAllele _)))
          .AD(alleles.map(allelicDepths.getOrElse(_, (0, 0))).map(_._1).toArray)
          .attribute("GGQ", item.germlineGenotypeQuality)
          .attribute("NPE", item.normalPercentError)
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
      }
    })

    val trigger = (
      (if (item.pooledCall) Seq("POOLED") else Seq.empty) ++
      (if (item.individualCall) Seq("INDIVIDUAL") else Seq.empty)).mkString("+")

    new VariantContextBuilder()
      .chr(item.contig)
      .start(adjustedPosition + 1) // one based based inclusive
      .stop(adjustedPosition + 1 + math.max(adjustAllele(item.ref).length - 1, 0))
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

    var prevChr = ""
    var prevStart = -1
    calls.foreach(call => {
      val variantContext = makeHtsjdVariantContext(call, sampleNames, reference)
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