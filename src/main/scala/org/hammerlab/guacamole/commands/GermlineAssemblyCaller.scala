package org.hammerlab.guacamole.commands

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.hammerlab.guacamole.Common.Arguments.GermlineCallerArgs
import org.hammerlab.guacamole.DistributedUtil.PerSample
import org.hammerlab.guacamole._
import org.hammerlab.guacamole.alignment.AlignmentState.AlignmentState
import org.hammerlab.guacamole.alignment.{AffineGapPenaltyAlignment, AlignmentState, ReadAlignment}
import org.hammerlab.guacamole.assembly.DeBruijnGraph
import org.hammerlab.guacamole.filters.GenotypeFilter.GenotypeFilterArguments
import org.hammerlab.guacamole.filters._
import org.hammerlab.guacamole.reads.{MappedRead, Read}
import org.hammerlab.guacamole.reference.{ReferenceBroadcast, ReferenceGenome}
import org.hammerlab.guacamole.variants.{Allele, AlleleConversions, AlleleEvidence, CalledAllele}
import org.hammerlab.guacamole.windowing.SlidingWindow
import org.kohsuke.args4j.{Option => Args4jOption}

/**
 *
 * Simple assembly based germline variant caller
 *
 */
object GermlineAssemblyCaller {

  class Arguments extends GermlineCallerArgs with GenotypeFilterArguments {

    @Args4jOption(name = "--kmer-size", usage = "Length of kmer used for DeBrujin Graph assembly")
    var kmerSize: Int = 45

    @Args4jOption(name = "--snv-window-range", usage = "Number of bases before and after to check for additional matches or deletions")
    var snvWindowRange: Int = 20

    @Args4jOption(name = "--min-average-base-quality", usage = "Minimum average of base qualities in the read")
    var minAverageBaseQuality: Int = 20

    @Args4jOption(name = "--min-alignment-quality", usage = "Minimum alignment qualities of the read")
    var minAlignmentQuality: Int = 30

    @Args4jOption(name = "--reference-fasta", required = false, usage = "Local path to a reference FASTA file")
    var referenceFastaPath: String = null

  }

  object Caller extends SparkCommand[Arguments] {
    override val name = "germline-assembly"
    override val description = "call germline variants by assembling the surrounding region of reads"

    /**
     * Extract alignments that are mismatches, insertions or deletions
     * @param alignment
     * @return Filtered set of alignments that are mismatches, insertions or deletions
     */
    def getVariantAlignments(alignment: ReadAlignment): Seq[(AlignmentState, Int)] = {
      alignment
        .alignments
        .zipWithIndex
        .filter(tup => AlignmentState.isGapAlignment(tup._1) || tup._1 == AlignmentState.Mismatch)
    }

    /**
     *
     * @param graph
     * @param windows
     * @param kmerSize
     * @param minAlignmentQuality
     * @param minAverageBaseQuality
     * @param minOccurrence
     * @param expectedPloidy
     * @return
     */
    def discoverHaplotypes(graph: Option[DeBruijnGraph],
                           windows: PerSample[SlidingWindow[MappedRead]],
                           kmerSize: Int,
                           minAlignmentQuality: Int,
                           minAverageBaseQuality: Int,
                           reference: ReferenceGenome,
                           minOccurrence: Int = 2,
                           expectedPloidy: Int = 2): (Option[DeBruijnGraph], Iterator[CalledAllele]) = {

      val currentWindow = windows(0)
      val locus = currentWindow.currentLocus

      val newReads = currentWindow.currentRegions()
      val filteredReads = newReads.filter(_.alignmentQuality > minAlignmentQuality)

      // TODO: Should update graph instead of rebuilding it
      // Need to keep track of reads removed from last update and reads added
      val currentGraph: DeBruijnGraph = DeBruijnGraph(
        filteredReads.map(_.sequence),
        kmerSize,
        minOccurrence
      )

      val referenceStart = (currentWindow.currentLocus - currentWindow.halfWindowSize).toInt
      val referenceEnd = (currentWindow.currentLocus + currentWindow.halfWindowSize).toInt

      val currentReference: Array[Byte] = reference.getReferenceSequence(
        currentWindow.referenceName,
        referenceStart,
        referenceEnd
      )

      def buildVariant(locus: Long, tuple: (AlignmentState, Int)): CalledAllele = {

        //val variantType = tuple._1
        val offset = tuple._2
        val allele = Allele(Seq(currentReference((locus + offset).toInt)), Bases.stringToBases("<ALT>"))
        val depth = filteredReads.length
        CalledAllele(
          currentWindow.newRegions.head.sampleName,
          currentWindow.newRegions.head.referenceContig,
          locus,
          allele,
          AlleleEvidence(
            likelihood = 1,
            readDepth = depth,
            alleleReadDepth = depth,
            forwardDepth = depth,
            alleleForwardDepth = depth,
            meanMappingQuality = 60,
            medianMappingQuality = 60,
            meanBaseQuality = 60,
            medianBaseQuality = 60,
            medianMismatchesPerRead = 0
          )
        )
      }

      val referenceKmerSource = currentReference.take(kmerSize)
      val referenceKmerSink = currentReference.takeRight(kmerSize)

      val paths = currentGraph.depthFirstSearch(
        referenceKmerSource,
        referenceKmerSink
      )

      // Take `ploidy` paths
      val topPaths = paths.take(expectedPloidy).map(DeBruijnGraph.mergeKmers)
      // Align paths to reference
      val alignments = topPaths.map(AffineGapPenaltyAlignment.align(_, currentReference))
      // Output variants
      val variantAlignments =
        alignments
          .flatMap(
            getVariantAlignments(_)
              .map(tup => buildVariant(locus, tup)
              )
          )
      (graph, variantAlignments.iterator)
    }

    override def run(args: Arguments, sc: SparkContext): Unit = {

      val readSet = Common.loadReadsFromArguments(
        args,
        sc,
        Read.InputFilters(mapped = true, nonDuplicate = true))

      Common.progress(
        "Loaded %,d mapped non-duplicate reads into %,d partitions.".format(readSet.mappedReads.count, readSet.mappedReads.partitions.length))

      val loci = Common.loci(args, readSet)
      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(args, loci, readSet.mappedReads)

      val kmerSize = args.kmerSize
      val minAlignmentQuality = args.minAlignmentQuality
      val minAverageBaseQuality = args.minAverageBaseQuality

      val lociOfInterest = DistributedUtil.pileupFlatMap[VariantLocus](
        readSet.mappedReads,
        lociPartitions,
        skipEmpty = true,
        function = pileup => VariantLocus(pileup).iterator,
        reference = None
      )

      val lociOfInterestSet = LociSet.union(
        lociOfInterest.map(
          locus => LociSet(locus.contig, locus.locus, locus.locus + 1))
          .collect(): _*
      )

      val lociOfInterestPartitions = DistributedUtil.partitionLociUniformly(
        args.parallelism,
        lociOfInterestSet
      )

      val reference = ReferenceBroadcast(args.referenceFastaPath, sc)
      val genotypes: RDD[CalledAllele] =
        DistributedUtil.windowFlatMapWithState[MappedRead, CalledAllele, Option[DeBruijnGraph]](
          Seq(readSet.mappedReads),
          lociPartitions,
          skipEmpty = true,
          halfWindowSize = args.snvWindowRange,
          initialState = None,
          (graph, window) => {
            discoverHaplotypes(
              graph,
              window,
              kmerSize,
              minAlignmentQuality,
              minAverageBaseQuality,
              reference)
          }
        )

      readSet.mappedReads.unpersist()

      val filteredGenotypes =
        GenotypeFilter(genotypes, args)
          .flatMap(AlleleConversions.calledAlleleToADAMGenotype)
      Common.writeVariantsFromArguments(args, filteredGenotypes)
      DelayedMessages.default.print()
    }

  }
}
