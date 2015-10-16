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
import org.hammerlab.guacamole.Bases
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
    val file = File.createTempFile("test-somatic-joint-caller", suffix)
    //file.deleteOnExit()
    file.getAbsolutePath
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
    val records = vcfRecords(reader)
    println("Guacamole calls: %,d. Gold calls: %,d.".format(records.length, recordsGold.length))



  }
}
