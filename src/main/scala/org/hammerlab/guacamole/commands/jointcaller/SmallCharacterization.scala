package org.hammerlab.guacamole.commands.jointcaller

import htsjdk.variant.variantcontext.VariantContext
import org.bdgenomics.adam.util.PhredUtils
import org.hammerlab.guacamole._
import org.hammerlab.guacamole.pileup.PileupElement
import org.hammerlab.guacamole.reference.ReferenceBroadcast

trait SmallCharacterization {
  val contig: String
  val position: Long
}

object SmallCharacterization {
  case class PileupStats(ref: String, elements: Seq[PileupElement]) {
    val alleleToElements = elements.sortBy(_.qualityScore * -1).groupBy(element => {
      Bases.basesToString(element.sequencedBases).toUpperCase
    })
    val alleleStats = alleleToElements
      .mapValues(items => {
        val positive = items.filter(_.read.isPositiveStrand)
        val negative = items.filter(!_.read.isPositiveStrand)
        val numToTake = Math.max(Math.min(positive.size, negative.size), 3)
        val adjustedItems = positive.take(numToTake) ++ negative.take(numToTake)
        (items.size, positive.size, adjustedItems)
      })
    val allelicDepths = alleleStats.mapValues(item => (item._1, item._2))
    val downsampledElements = alleleStats.mapValues(_._3).values.flatten

    val nonRefAlleles = allelicDepths.filterKeys(_ != ref).toSeq.sortBy(_._2._1 * -1).map(_._1)

    val topAlt = nonRefAlleles.headOption.getOrElse("N")
    val secondAlt = if (nonRefAlleles.size > 1) nonRefAlleles(1) else "N"

    lazy val readsSupportingAllele = alleleToElements
      .mapValues(_.filter(_.read.alignmentQuality > 0).map(_.read.name).toSet)
      .withDefaultValue(Set.empty)

    def logLikelihoodPileup(mixture: Map[String, Double]): Double = {
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
      downsampledElements.map(logLikelihoodPileupElement _).sum
    }
  }
}
