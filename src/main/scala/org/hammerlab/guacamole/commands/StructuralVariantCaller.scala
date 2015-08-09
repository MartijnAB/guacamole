package org.hammerlab.guacamole.commands

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.hammerlab.guacamole.{ DistributedUtil, SparkCommand, Common }
import org.hammerlab.guacamole.reads.{ PairedMappedRead, Read }
import org.kohsuke.args4j.{ Option => Args4JOption }

import scalax.collection.Graph
import scalax.collection.edge._
import scalax.collection.mutable.{ Graph => MutableGraph }
import scalax.collection.GraphPredef._
import scalax.collection.edge.Implicits._

/**
 * Structural Variant caller
 *
 * This currently looks for regions with abnormally large insert sizes. This typically indicates
 * that there's been a large deletion.
 *
 * It does this by finding the average insert size of all mate pairs which overlap each 25-base pair "block" of the
 * genome. If this average is significantly different than the average insert size genome-wide, then it's evidence that
 * there may be a structural variant.
 *
 * The output is a text file of genomic ranges where there's some evidence of an SV.
 */
object StructuralVariant {

  val MAX_INSERT_SIZE = 25000
  val BLOCK_SIZE = 25

  protected class Arguments extends Common.Arguments.Reads with DistributedUtil.Arguments {
    @Args4JOption(name = "--filter-contig", usage = "Filter to alignments where either mate is in this contig.")
    var filterContig: String = ""

    @Args4JOption(name = "--output", usage = "Path to output text file.")
    var output: String = ""
  }

  case class Location(
    contig: String,
    position: Long)

  case class GenomeRange(
    contig: String,
    start: Long,
    stop: Long)

  case class MedianStats(
    median: Double,
    mad: Double)

  // An undirected, weighted graph of evidence for structural variants.
  // Nodes are paired reads.
  // Edges represent compatibility between the induced structural variants,
  // with the weight being the strength of compatibility.
  type PairedReadGraph = Graph[PairedMappedRead, WUnDiEdge]

  object Caller extends SparkCommand[Arguments] {
    override val name = "structural-variant"
    override val description = "Find structural variants, e.g. large deletions."

    // The sign of inferredInsertSize follows the orientation of the read. This returns
    // an insert size which is positive in the common case, regardless of the orientation.
    def orientedInsertSize(r: PairedMappedRead): Int = {
      val sgn = if (r.read.isPositiveStrand) { 1 } else { -1 }
      r.inferredInsertSize * sgn
    }

    // Extract median and MAD (Median Absolute Deviation) from an unordered Array.
    def medianStats[T <% Double](xs: Array[T])(implicit ord: Ordering[T]): MedianStats = {

      // Extract the median of a sorted, non-empty array of numbers
      def getMedian[T <% Double](nums: Array[T]): Double = {
        val len = nums.length
        if (len % 2 == 0) {
          0.5 * (nums(len / 2 - 1) + nums(len / 2))
        } else {
          1.0 * nums(len / 2)
        }
      }

      if (xs.isEmpty) {
        MedianStats(0.0, 0.0)
      } else {
        val nums = xs.sorted
        val median = getMedian(nums)
        val residuals = nums.map(x => Math.abs(1.0 * x - median)).sorted
        val mad = getMedian(residuals)

        MedianStats(median, mad)
      }
    }

    case class ExceptionalReadsReturnType(
      readsInRange: RDD[PairedMappedRead],
      insertSizes: RDD[Int],
      insertStats: MedianStats,
      maxNormalInsertSize: Int,
      exceptionalReads: RDD[PairedMappedRead])

    def getExceptionalReads(reads: RDD[PairedMappedRead]): ExceptionalReadsReturnType = {
      // Pare down to mate pairs which aren't highly displaced & aren't inverted or translocated.
      val readsInRange = reads.filter(r => {
        val s = r.inferredInsertSize
        r.read.referenceContig == r.mate.referenceContig &&
          r.read.isPositiveStrand != r.mate.isPositiveStrand &&
          s < MAX_INSERT_SIZE
      })
      readsInRange.persist()

      val insertSizes = readsInRange.map(orientedInsertSize)
      val insertStats = medianStats(insertSizes.takeSample(false, 100000))
      println("insert stats: " + insertStats)

      // Find blocks where the average insert size differs significantly from the genome-wide median.
      // Following DELLY, we use a threshold of 5 median absolute deviations above the median.
      val maxNormalInsertSize = (insertStats.median + 5 * insertStats.mad).toInt

      val exceptionalReads = readsInRange.filter(_.inferredInsertSize > maxNormalInsertSize)

      ExceptionalReadsReturnType(
        readsInRange,
        insertSizes,
        insertStats,
        maxNormalInsertSize,
        exceptionalReads
      )
    }

    // Given two pairs of reads, is there a deletion which would make both of them have normal insert sizes?
    def areReadsCompatible(read1: PairedMappedRead, read2: PairedMappedRead, maxNormalInsertSize: Long): Boolean = {
      if (read1.minPos > read2.minPos) {
        areReadsCompatible(read2, read1, maxNormalInsertSize)
      } else {
        // This matches the DELLY logic; the last five lines are copy/paste. See
        // https://github.com/tobiasrausch/delly/blob/7464cacfe1d420ac5bd7ed40a3093a80a93a6964/src/tags.h#L441-L445
        val (pair1Min, pair1GapMin, pair1GapMax, pair1Max) = read1.startsAndStops
        val pair1ReadLength = read1.readLength

        val (pair2Min, pair2GapMin, pair2GapMax, pair2Max) = read2.startsAndStops
        val pair2ReadLength = read2.readLength

        !(((pair2GapMin - pair1Min) > maxNormalInsertSize) ||
          ((pair2GapMax < pair1GapMax) && ((pair1Max - pair2GapMax) > maxNormalInsertSize)) ||
          ((pair2GapMax >= pair1GapMax) && ((pair2Max - pair1GapMax) > maxNormalInsertSize)) ||
          ((pair1GapMax < pair2Min) || (pair2GapMax < pair1Min)))
      }
    }

    // Given a list of reads with abnormally large inserts, construct a graph representation of the structural variants
    // where:
    // Vertices = paired reads
    // Edges = two sets of paired reads could represent the same structural variant
    // For more details, see the DELLY paper.
    def buildVariantGraph(exceptionalReads: Iterable[PairedMappedRead], maxNormalInsertSize: Int): PairedReadGraph = {
      val reads = exceptionalReads.toArray.sortBy(_.minPos)
      val graph = MutableGraph[PairedMappedRead, WUnDiEdge]()

      for { (read, i) <- reads.zipWithIndex } {
        val start = read.minPos
        var j = i + 1
        var done = false
        // Loop through subsequent reads until it's impossible to find another compatible one.
        while (j < reads.length && !done) {
          val nextRead = reads(j)
          val nextStart = nextRead.minPos
          if (Math.abs(nextStart + nextRead.readLength - start) > maxNormalInsertSize) {
            done = true
          } else {
            if (areReadsCompatible(read, nextRead, maxNormalInsertSize)) {
              graph += (read ~ nextRead) % 1
            }
          }
          j += 1
        }
      }

      // TODO: exclude isolated reads from the graph
      graph
    }

    override def run(args: Arguments, sc: SparkContext): Unit = {
      val readSet = Common.loadReadsFromArguments(args, sc, Read.InputFilters(nonDuplicate = true))
      val pairedMappedReads = readSet.mappedPairedReads
      val firstInPair = pairedMappedReads.filter(_.isFirstInPair).flatMap(PairedMappedRead(_))

      val ExceptionalReadsReturnType(_, _, _, maxNormalInsertSize, exceptionalReads) = getExceptionalReads(firstInPair)

      val readsByContig = exceptionalReads.groupBy(_.read.referenceContig)
      val svGraph = readsByContig.mapValues(buildVariantGraph(_, maxNormalInsertSize))

      svGraph.coalesce(1).saveAsTextFile(args.output)
    }
  }
}