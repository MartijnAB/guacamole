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

package org.hammerlab.guacamole.commands.jointcaller

import java.util

import htsjdk.samtools.util.Locus
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.variantcontext.{ Allele, GenotypeBuilder, VariantContext, VariantContextBuilder }
import htsjdk.variant.vcf._
import org.apache.http.client.utils.URLEncodedUtils
import org.apache.spark.SparkContext
import org.bdgenomics.adam.util.PhredUtils
import org.hammerlab.guacamole.Common.Arguments.NoSequenceDictionary
import org.hammerlab.guacamole._
import org.hammerlab.guacamole.pileup._
import org.hammerlab.guacamole.reads._
import org.hammerlab.guacamole.reference.ReferenceBroadcast
import org.kohsuke.args4j.{ Argument, Option => Args4jOption }

import scala.collection.mutable.ArrayBuffer
import scala.collection.{ JavaConversions, mutable }

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

  object Caller extends SparkCommand[Arguments] {
    override val name = "somatic-joint"
    override val description = "somatic caller for any number of samples from the same patient"

    override def run(args: Arguments, sc: SparkContext): Unit = {
      val inputs = Inputs.parseMultiple(args.inputs)

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
      val groupedInputs = Inputs.GroupedInputs(inputs)

      // TODO: we currently are re-shuffling on every pileupFlatMap call.
      // Call germline small variants.
      val germlineCalls = DistributedUtil.pileupFlatMapMultipleRDDs(
        groupedInputs.normalDNA.map(readSets(_).mappedReads),
        lociPartitions,
        false, // skip empty
        pileups => {
          val characterization = GenomicCharacterization.GermlineCharacterization(groupedInputs.normalDNA.map(pileups(_)))
          val emit = (characterization.isVariant
            || broadcastForceCallLoci.value.onContig(pileups(0).referenceName).contains(pileups(0).locus))
          if (emit) Iterator(characterization) else Iterator.empty
        }, referenceGenome = reference).collect

      Common.progress("Called %,d germline variants.".format(germlineCalls.length))

      if (args.outSmallGermlineVariants.nonEmpty) {
        Common.progress("Writing germline variants.")
        GenomicCharacterization.GermlineCharacterization.writeVcf(
          args.outSmallGermlineVariants, args.normalSampleName, readSets, germlineCalls, reference.get)
        Common.progress("Wrote %,d calls to %s".format(germlineCalls.length, args.outSmallGermlineVariants))
      }

      /*
      val emitGermlineCalls = args.outSmallGermlineVariants.nonEmpty
      val calls = DistributedUtil.pileupFlatMapMultipleRDDs(
        readSets.map(_.mappedReads),
        lociPartitions,
        false, // skip empty
        pileups => {
          val germlineCharacterization = GenomicCharacterization.GermlineCharacterization(groupedInputs.normalDNA.map(pileups(_)))
          val emitGermlineCall = (
            emitGermlineCalls && (
              germlineCharacterization.isVariant ||
              broadcastForceCallLoci.value.onContig(pileups(0).referenceName).contains(pileups(0).locus)))
          if (emitGermlineCall) Iterator(germlineCharacterization) else Iterator.empty
        }, referenceGenome = reference).collect
       */

    }

  }

  def combinedPileup(pileups: Seq[Pileup]) = {
    val elements = pileups.flatMap(_.elements)
    Pileup(pileups(0).referenceName, pileups(0).locus, pileups(0).referenceBase, elements)
  }
}
