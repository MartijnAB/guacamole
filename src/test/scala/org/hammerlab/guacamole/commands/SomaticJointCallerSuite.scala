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

import au.com.bytecode.opencsv.CSVParser
import htsjdk.variant.variantcontext.{Allele, VariantContextBuilder, VariantContext}
import htsjdk.variant.vcf.VCFFileReader
import org.hammerlab.guacamole.commands.jointcaller.SomaticJoint
import SomaticJoint.Arguments
import org.hammerlab.guacamole.util.{ TestUtil, GuacFunSuite }
import org.hammerlab.guacamole.{ LociMap, Bases }
import org.hammerlab.guacamole.filters.SomaticGenotypeFilter
import org.hammerlab.guacamole.pileup.Pileup
import org.hammerlab.guacamole.reads.MappedRead
import org.scalatest.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import scala.collection.JavaConversions
import scala.collection.mutable.ArrayBuffer
import scala.io.Source

class SomaticJointCallerSuite extends GuacFunSuite with Matchers {

  def loadPileup(filename: String, referenceName: String, locus: Long = 0): Pileup = {
    val records = TestUtil.loadReads(sc, filename).mappedReads
    val localReads = records.collect
    Pileup(localReads, referenceName, locus)
  }

  var tempFileNum = 0
  def tempFile(suffix: String): String = {
    tempFileNum += 1
    "/tmp/test-somatic-joint-caller-%d.vcf".format(tempFileNum)
    //val file = File.createTempFile("test-somatic-joint-caller", suffix)
    //file.deleteOnExit()
    //file.getAbsolutePath
  }

  val na12878SubsetBam = TestUtil.testDataPath(
    "illumina-platinum-na12878-extract/NA12878.10k_variants.plus_chr1_3M-3.1M.chr_fixed.bam")

  val na12878ExpectedCallsVCF = TestUtil.testDataPath(
    "illumina-platinum-na12878-extract/NA12878.subset.vcf")

  val chr1PrefixFasta = TestUtil.testDataPath("illumina-platinum-na12878-extract/chr1.prefix.fa")

  val cancerWGS1Bams = Seq("normal.bam", "primary.bam", "recurrence.bam").map(
    name => TestUtil.testDataPath("cancer-wgs1/" + name))

  val cancerWGS1ExpectedSomaticCallsCSV = TestUtil.testDataPath("cancer-wgs1/variants.csv")

  def vcfRecords(reader: VCFFileReader): Seq[VariantContext] = {
    val results = new ArrayBuffer[VariantContext]
    val iterator = reader.iterator
    while (iterator.hasNext) {
      val value = iterator.next()
      results += value
    }
    results
  }

  case class VariantFromVarlensCSV(
                                    genome: String,
                                    contig: String,
                                    interbaseStart: Int,
                                    interbaseEnd: Int,
                                    ref: String,
                                    alt: String,
                                    tumor: String,
                                    normal: String,
                                    validation: String)
  {}

  def csvRecords(filename: String): Seq[VariantFromVarlensCSV] = {
    Source.fromFile(filename).getLines.map(_.split(",").map(_.trim).toSeq).toList match {
      case header :: records => {
        records.map(record => {
          val fields = header.zip(record).toMap
          VariantFromVarlensCSV(
            fields("genome"),
            fields("contig"),
            fields("interbase_start").toInt,
            fields("interbase_end").toInt,
            fields("ref"),
            fields("alt"),
            fields("tumor"),
            fields("normal"),
            fields("validation"))
        })
      }
    }
  }

  case class VCFComparison(expected: Seq[VariantContext], experimental: Seq[VariantContext]) {

    val mapExpected = VCFComparison.makeLociMap(expected)
    val mapExperimental = VCFComparison.makeLociMap(experimental)

    val exactMatch = new ArrayBuffer[(VariantContext, VariantContext)]
    val partialMatch = new ArrayBuffer[(VariantContext, VariantContext)]
    val uniqueToExpected = new ArrayBuffer[VariantContext]
    val uniqueToExperimental = new ArrayBuffer[VariantContext]

    VCFComparison.accumulate(expected, mapExperimental, exactMatch, partialMatch, uniqueToExpected)

    {
      // Accumulate result.{exact,partial}Match in throwaway arrays so we don't double count.
      val exactMatch2 = new ArrayBuffer[(VariantContext, VariantContext)]
      val partialMatch2 = new ArrayBuffer[(VariantContext, VariantContext)]

      VCFComparison.accumulate(experimental, mapExpected, exactMatch2, partialMatch2, uniqueToExperimental)
      // assert(exactMatch2.size == exactMatch.size)
      // assert(partialMatch2.size == partialMatch.size)
    }

    def sensitivity = exactMatch.size * 100.0 / expected.size
    def specificity = exactMatch.size * 100.0 / experimental.size

    def summary(): String = {
      /*
      def show(count: Int): String = {
        "%,d (%0.2f%%)".format()
      }
      */

      Seq(
        "exact match: %,d".format(exactMatch.size),
        "partial match: %,d".format(partialMatch.size),
        "unique to expected: %,d".format(uniqueToExpected.size),
        "unique to experimental: %,d".format(uniqueToExperimental.size),
        "sensitivity: %1.2f%%".format(sensitivity),
        "specificity: %1.2f%%".format(specificity)
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
        map.onContig(record1.getContig).get(record1.getStart) match {
          case Some(record2) => {
            if (variantToString(record1) == variantToString(record2)) {
              exactMatch += ((record1, record2))
            } else {
              partialMatch += ((record1, record2))
            }
          }
          case None => unique += record1
        }
      })
    }

    private def makeLociMap(records: Seq[VariantContext]): LociMap[VariantContext] = {
      val builder = LociMap.newBuilder[VariantContext]()
      records.foreach(record => {
        // Switch from zero based inclusive to interbase coordinates.
        builder.put(record.getContig, record.getStart, record.getEnd + 1, record)
      })
      builder.result
    }
  }

  def variantToString(variant: VariantContext, verbose: Boolean = false): String = {
    val genotype = variant.getGenotype(0)

    if (verbose) {
      variant.toString
    } else {
      "%s:%d-%d %s > %s %s".format(
        variant.getContig,
        variant.getStart,
        variant.getEnd,
        variant.getReference,
        JavaConversions.collectionAsScalaIterable(variant.getAlternateAlleles).map(_.toString).mkString(","),
        genotype.getType.toString)
    }
  }

  def printSamplePairs(pairs: Seq[(VariantContext, VariantContext)], num: Int = 20): Unit = {
    val sample = pairs.take(num)
    println("Showing %,d / %,d.".format(sample.size, pairs.size))
    sample.zipWithIndex.foreach({
      case (pair, num) => {
        println("(%4d) %20s vs %20s \tDETAILS: %20s vs %20s".format(
          num + 1,
          variantToString(pair._1, false),
          variantToString(pair._2, false),
          variantToString(pair._1, true),
          variantToString(pair._2, true)))
      }
    })
  }

  def printSample(items: Seq[VariantContext], num: Int = 20): Unit = {
    val sample = items.take(num)
    println("Showing %,d / %,d.".format(sample.size, items.size))
    sample.zipWithIndex.foreach({
      case (item, num) => {
        println("(%4d) %20s \tDETAILS: %29s".format(
          num + 1,
          variantToString(item, false),
          variantToString(item, true)))
      }
    })
  }

  def compareToCSV(experimentalFile: String, expectedFile: String): Unit = {
    val recordsExpected = csvRecords(expectedFile)
    val readerExperimental = new VCFFileReader(new File(experimentalFile), false)
    val recordsExperimental = vcfRecords(readerExperimental)

    println("Experimental calls: %,d. Gold calls: %,d.".format(recordsExperimental.length, recordsExpected.length))

  }

  def compareToVCF(experimentalFile: String, expectedFile: String): Unit = {

    val readerExpected = new VCFFileReader(new File(expectedFile), false)
    val recordsExpected = vcfRecords(readerExpected)
    val reader = new VCFFileReader(new File(experimentalFile), false)
    val recordsGuacamole = vcfRecords(reader)

    println("Experimental calls: %,d. Gold calls: %,d.".format(recordsGuacamole.length, recordsExpected.length))

    def onlyIndels(calls: Seq[VariantContext]): Seq[VariantContext] = {
      calls.filter(call => call.getReference.length != 1 ||
        JavaConversions.collectionAsScalaIterable(call.getAlternateAlleles).exists(_.length != 1))
    }

    val comparisonFull = VCFComparison(recordsExpected, recordsGuacamole)
    val comparisonPlatinumOnly = VCFComparison(
      recordsExpected.filter(_.getAttributeAsString("metal", "") == "platinum"),
      recordsGuacamole)
    val comparisonFullIndels = VCFComparison(onlyIndels(recordsExpected), onlyIndels(recordsGuacamole))

    println("Sensitivity on platinum: %f".format(comparisonPlatinumOnly.sensitivity))
    println("Specificity on full: %f".format(comparisonFull.specificity))

    println(comparisonFull.summary)

    println("MISSED CALLS IN PLATINUM WITH DEPTH")
    printSamplePairs(comparisonPlatinumOnly.partialMatch.filter(
      pair => pair._2.getGenotype(0).isHomRef && pair._2.getGenotype(0).getDP > 5))
    println()

    println("BAD CALLS")
    printSamplePairs(comparisonFull.partialMatch.filter(
      pair => !pair._2.getGenotype(0).isHomRef))
    println()

    println("INDEL PERFORMANCE")
    println(comparisonFullIndels.summary)

    println("INDEL EXACT MATCHES")
    printSamplePairs(comparisonFullIndels.exactMatch)

    println("INDEL MISSED CALLS WITH DEPTH")
    printSamplePairs(comparisonFullIndels.partialMatch.filter(
      pair => pair._2.getGenotype(0).isHomRef && pair._2.getGenotype(0).getDP > 5))
    println()

    println("INDEL BAD CALLS")
    printSamplePairs(comparisonFullIndels.partialMatch.filter(
      pair => !pair._2.getGenotype(0).isHomRef))
    println()

  }

  sparkTest("somatic calling on subset of 3-sample cancer patient 1") {
    val germlineResultFile = tempFile(".vcf")
    val somaticResultFile = tempFile(".vcf")
    println(germlineResultFile, somaticResultFile)

    if (true) {
      val args = new SomaticJoint.Arguments()
      args.outSmallGermlineVariants = germlineResultFile
      args.outSmallSomaticVariants = somaticResultFile
      args.referenceFastaPath = "/Users/tim/sinai/data/ucsc.hg19.fasta.gz"
      args.loci = ((1).until(22).map(i => "chr%d".format(i)) ++ Seq("chrX", "chrY")).mkString(",")
      args.loci = "chr4:71115210-71115219"

      args.inputs = cancerWGS1Bams.toArray
      //args.forceCallLociFromFile = na12878ExpectedCallsVCF
      SomaticJoint.Caller.run(args, sc)
    }

    println("************* CANCER WGS1 SOMATIC CALLS *************")

    compareToCSV(somaticResultFile, cancerWGS1ExpectedSomaticCallsCSV)


  }

  sparkTest("germline calling on subset of illumina platinum NA12878") {
    if (false) {
      val resultFile = tempFile(".vcf")
      println(resultFile)

      if (true) {
        val args = new SomaticJoint.Arguments()
        args.outSmallGermlineVariants = resultFile
        args.inputs = Seq(na12878SubsetBam).toArray
        args.loci = "chr1:0-6700000"
        args.forceCallLociFromFile = na12878ExpectedCallsVCF
        args.referenceFastaPath = chr1PrefixFasta
        SomaticJoint.Caller.run(args, sc)
      }

      println("************* GUACAMOLE *************")
      compareToVCF(resultFile, na12878ExpectedCallsVCF)

      println("************* UNIFIED GENOTYPER *************")
      compareToVCF(TestUtil.testDataPath(
        "illumina-platinum-na12878-extract/unified_genotyper.vcf"),
        na12878ExpectedCallsVCF)

      println("************* HaplotypeCaller *************")
      compareToVCF(TestUtil.testDataPath(
        "illumina-platinum-na12878-extract/haplotype_caller.vcf"),
        na12878ExpectedCallsVCF)
    }
  }
}
