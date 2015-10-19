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

import java.io.File

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import org.hammerlab.guacamole.commands.SomaticJoint.Arguments
import org.hammerlab.guacamole.util.{ TestUtil, GuacFunSuite }
import org.hammerlab.guacamole.{LociMap, Bases}
import org.hammerlab.guacamole.filters.SomaticGenotypeFilter
import org.hammerlab.guacamole.pileup.Pileup
import org.hammerlab.guacamole.reads.MappedRead
import org.scalatest.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import scala.collection.JavaConversions
import scala.collection.mutable.ArrayBuffer

class SomaticJointCallerSuite extends GuacFunSuite with Matchers with TableDrivenPropertyChecks {

  def loadPileup(filename: String, referenceName: String, locus: Long = 0): Pileup = {
    val records = TestUtil.loadReads(sc, filename).mappedReads
    val localReads = records.collect
    Pileup(localReads, referenceName, locus)
  }

  def tempFile(suffix: String): String = {
    "/tmp/test-somatic-joint-caller.vcf"
    //val file = File.createTempFile("test-somatic-joint-caller", suffix)
    //file.deleteOnExit()
    //file.getAbsolutePath
  }

  val na12878_subset_bam = TestUtil.testDataPath(
    "illumina-platinum-na12878-extract/NA12878.10k_variants.plus_chr1_3M-3.1M.bam")

  val na12878_gold_calls_vcf = TestUtil.testDataPath(
    "illumina-platinum-na12878-extract/NA12878.subset.vcf")

  def vcfRecords(reader: VCFFileReader): Seq[VariantContext] = {
    val results = new ArrayBuffer[VariantContext]
    val iterator = reader.iterator
    while (iterator.hasNext) {
      val value = iterator.next()
      results += value
    }
    results
  }


  case class VCFComparison(gold: Seq[VariantContext], experimental: Seq[VariantContext]) {

    val mapGold = VCFComparison.makeLociMap(gold)
    val mapExperimental = VCFComparison.makeLociMap(experimental)

    val exactMatch = new ArrayBuffer[(VariantContext, VariantContext)]
    val partialMatch = new ArrayBuffer[(VariantContext, VariantContext)]
    val uniqueToGold = new ArrayBuffer[VariantContext]
    val uniqueToExperimental = new ArrayBuffer[VariantContext]

    VCFComparison.accumulate(gold, mapExperimental, exactMatch, partialMatch, uniqueToGold)

    {
      // Accumulate result.{exact,partial}Match in throwaway arrays so we don't double count.
      val exactMatch2 = new ArrayBuffer[(VariantContext, VariantContext)]
      val partialMatch2 = new ArrayBuffer[(VariantContext, VariantContext)]

      VCFComparison.accumulate(experimental, mapGold, exactMatch2, partialMatch2, uniqueToExperimental)
      // assert(exactMatch2.size == exactMatch.size)
      // assert(partialMatch2.size == partialMatch.size)
    }

    def summary(): String = {
      /*
      def show(count: Int): String = {
        "%,d (%0.2f%%)".format()
      }
      */

      Seq(
        "exact match: %,d".format(exactMatch.size),
        "partial match: %,d".format(partialMatch.size),
        "unique to gold: %,d".format(uniqueToGold.size),
        "unique to experimental: %,d".format(uniqueToExperimental.size),
        "sensitivity (exact): %1.2f%%".format(exactMatch.size * 100.0 / gold.size),
        "sensitivity (partial): %1.2f%%".format((exactMatch.size + partialMatch.size) * 100.0 / gold.size),
        "specificity (exact): %1.2f%%".format(exactMatch.size * 100.0 / experimental.size),
        "specificity (partial): %1.2f%%".format((exactMatch.size + partialMatch.size) * 100.0 / experimental.size)
      ).mkString("\n")
    }
    
  }
  object VCFComparison {
    private def accumulate(
                    records: Seq[VariantContext],
                    map: LociMap[VariantContext],
                    exactMatch: ArrayBuffer[(VariantContext, VariantContext)],
                    partialMatch: ArrayBuffer[(VariantContext, VariantContext)],
                    unique: ArrayBuffer[VariantContext]): Unit = {
      records.foreach(record1 => {
        map.onContig(record1.getChr).get(record1.getStart) match {
          case Some(record2) => {
            if (variantToString(record1) == variantToString(record2)) {
              exactMatch += ((record1, record2))
            } else {
              partialMatch += ((record1, record2))
            }
          }
          case None => (unique += record1)
        }
      })
    }

    private def makeLociMap(records: Seq[VariantContext]): LociMap[VariantContext] = {
      val builder = LociMap.newBuilder[VariantContext]()
      records.foreach(record => {
        // Switch from zero based inclusive to interbase coordinates.
        builder.put(record.getChr, record.getStart, record.getEnd + 1, record)
      })
      builder.result
    }
  }

  def variantToString(variant: VariantContext): String = {
    val genotype = variant.getGenotype(0)

    "%s:%d-%d %s > %s %s".format(
      variant.getChr,
      variant.getStart,
      variant.getEnd,
      variant.getReference,
      JavaConversions.collectionAsScalaIterable(variant.getAlternateAlleles).map(_.toString).mkString(","),
      genotype.getType.toString)
  }

  def printSamplePairs(pairs: Seq[(VariantContext, VariantContext)], num: Int = 20): Unit = {
    val sample = pairs.take(num)
    println("Showing %,d / %,d.".format(sample.size, pairs.size))
    sample.zipWithIndex.foreach({case (pair, num) => {
      println("(%4d) %20s vs %20s".format(num + 1, variantToString(pair._1), variantToString(pair._2)))
    }})
  }

  /*
  def printSample(items: Seq[VariantContext], num: Int = 20): Unit = {
    val sample = pairs.take(num)
    println("Showing %,d / %,d.".format(sample))
    sample.foreach(pair => {
      println("%20s vs %20s".format(variantToString(pair._1), variantToString(pair._2)))
    })
  }
  */

  sparkTest("germline calling on subset of illumina platinum NA12878") {
    val resultFile = tempFile(".vcf")
    println(resultFile)

    val args = new SomaticJoint.Arguments()
    args.outSmallGermlineVariants = resultFile
    args.inputs = Seq(na12878_subset_bam).toArray
    args.loci = "chr1:0-5000000"
    SomaticJoint.Caller.run(args, sc)

    val readerGold = new VCFFileReader(new File(na12878_gold_calls_vcf), false)
    val recordsGold = vcfRecords(readerGold)
    val reader = new VCFFileReader(new File(resultFile), false)
    val recordsGuacamole = vcfRecords(reader)
    println("Guacamole calls: %,d. Gold calls: %,d.".format(recordsGuacamole.length, recordsGold.length))
    
    val comparison = VCFComparison(recordsGold, recordsGuacamole)
    println(comparison.summary)

    println("EXACT MATCHES")
    printSamplePairs(comparison.exactMatch)
    println()

    println("PARTIAL MATCHES")
    printSamplePairs(comparison.partialMatch)
    println()


    println(comparison.summary)


  }
}
