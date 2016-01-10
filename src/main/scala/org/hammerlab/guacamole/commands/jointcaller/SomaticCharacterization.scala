package org.hammerlab.guacamole.commands.jointcaller

import htsjdk.variant.variantcontext.{VariantContextBuilder, GenotypeBuilder, Allele, VariantContext}
import org.hammerlab.guacamole.DistributedUtil._
import org.hammerlab.guacamole.ReadSet
import org.hammerlab.guacamole.commands.jointcaller.GenomicCharacterization.PileupStats
import org.hammerlab.guacamole.pileup.Pileup
import org.hammerlab.guacamole.reference.ReferenceBroadcast
import sun.reflect.generics.reflectiveObjects.NotImplementedException

import scala.collection.JavaConversions

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

    val pooledPosteriors = unnormalizedPosteriors(pooledStats)
    val perSamplePosteriors = sampleStats.map(unnormalizedPosteriors _)
    val noCallGenotype = Map(ref -> 1.0)

    SomaticCharacterization(
      tumorDNAPileups(0).referenceName,
      tumorDNAPileups(0).locus,
      ref,
      topAlt,
      call(pooledPosteriors) != noCallGenotype || perSamplePosteriors.exists(posterior => call(posterior) != noCallGenotype),
      sampleStats.map(_.allelicDepths))
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