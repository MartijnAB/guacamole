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

package org.hammerlab.guacamole.pileup

import org.hammerlab.guacamole.Bases
import org.hammerlab.guacamole.reads.MappedRead
import org.hammerlab.guacamole.util.TestUtil.Implicits._
import org.hammerlab.guacamole.util.TestUtil.assertBases
import org.hammerlab.guacamole.util.{ GuacFunSuite, TestUtil }
import org.hammerlab.guacamole.variants.Allele
import org.scalatest.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

class PileupSuite extends GuacFunSuite with Matchers with TableDrivenPropertyChecks {

  //lazy so that this is only accessed from inside a spark test where SparkContext has been initialized
  lazy val testAdamRecords = TestUtil.loadReads(sc, "different_start_reads.sam").mappedReads.collect()

  // The follow two functions provide a convenience for testing but not efficient as they rebuild
  // the reference sequence for a read every time
  def advancePileupElement(element: PileupElement, locus: Long): PileupElement = {
    element.advanceToLocus(locus, element.read.getReferenceBaseAtLocus(locus))
  }

  def pileupElementFromRead(read: MappedRead, locus: Long): PileupElement = {
    PileupElement(read, locus, read.getReferenceBaseAtLocus(locus))
  }

  def loadPileup(filename: String, locus: Long = 0): Pileup = {
    val records = TestUtil.loadReads(sc, filename).mappedReads
    val localReads = records.collect
    Pileup(localReads, localReads(0).referenceContig, locus)
  }

  sparkTest("create pileup from long insert reads") {
    val reads = Seq(
      TestUtil.makeRead("TCGATCGA", "8M", "8", 1),
      TestUtil.makeRead("TCGATCGA", "8M", "8", 1),
      TestUtil.makeRead("TCGACCCTCGA", "4M3I4M", "8", 1))

    val noPileup = Pileup(reads, "chr1", 0).elements
    assert(noPileup.size === 0)

    val firstPileup = Pileup(reads, "chr1", 1)
    firstPileup.elements.forall(_.isMatch) should be(true)
    firstPileup.elements.forall(_.qualityScore == 31) should be(true)

    val insertPileup = Pileup(reads, "chr1", 4)
    insertPileup.elements.exists(_.isInsertion) should be(true)
    insertPileup.elements.forall(_.qualityScore == 31) should be(true)

    insertPileup.elements(0).alignment should equal(Match('A', 31.toByte))
    insertPileup.elements(1).alignment should equal(Match('A', 31.toByte))
    insertPileup.elements(2).alignment should equal(Insertion("ACCC", Seq(31, 31, 31, 31).map(_.toByte)))
  }

  sparkTest("create pileup from long insert reads; different qualities in insertion") {
    val reads = Seq(
      TestUtil.makeRead("TCGATCGA", "8M", "8", 1, "chr1", Some(Seq(10, 15, 20, 25, 10, 15, 20, 25))),
      TestUtil.makeRead("TCGATCGA", "8M", "8", 1, "chr1", Some(Seq(10, 15, 20, 25, 10, 15, 20, 25))),
      TestUtil.makeRead("TCGACCCTCGA", "4M3I4M", "8", 1, "chr1", Some(Seq(10, 15, 20, 25, 5, 5, 5, 10, 15, 20, 25))))

    val insertPileup = Pileup(reads, "chr1", 4)
    insertPileup.elements.exists(_.isInsertion) should be(true)
    insertPileup.elements.exists(_.qualityScore == 5) should be(true)

    insertPileup.elements.forall(
      _.alignment match {
        case Match(_, quality)       => quality == 25
        case Insertion(_, qualities) => qualities.sameElements(Seq(25, 5, 5, 5))
        case _                       => false
      }) should be(true)
  }

  sparkTest("create pileup from long insert reads, right after insertion") {
    val reads = Seq(
      TestUtil.makeRead("TCGATCGA", "8M", "8", 1, "chr1", Some(Seq(10, 15, 20, 25, 10, 15, 20, 25))),
      TestUtil.makeRead("TCGATCGA", "8M", "8", 1, "chr1", Some(Seq(10, 15, 20, 25, 10, 15, 20, 25))),
      TestUtil.makeRead("TCGACCCTCGA", "4M3I4M", "8", 1, "chr1", Some(Seq(10, 15, 20, 25, 5, 5, 5, 10, 15, 20, 25)))
    )

    val noPileup = Pileup(reads, "chr1", 0).elements
    noPileup.size should be(0)

    val pastInsertPileup = Pileup(reads, "chr1", 5)
    pastInsertPileup.elements.foreach(_.isMatch should be(true))

    pastInsertPileup.elements.foreach(_.qualityScore should be(10))

  }

  sparkTest("create pileup from long insert reads; after insertion") {
    val reads = Seq(
      TestUtil.makeRead("TCGATCGA", "8M", "8", 1),
      TestUtil.makeRead("TCGATCGA", "8M", "8", 1),
      TestUtil.makeRead("TCGACCCTCGA", "4M3I4M", "8", 1))
    val lastPileup = Pileup(reads, "chr1", 7)
    lastPileup.elements.foreach(e => assertBases(e.sequencedBases, "G"))
    lastPileup.elements.forall(_.isMatch) should be(true)
  }

  sparkTest("create pileup from long insert reads; end of read") {

    val reads = Seq(
      TestUtil.makeRead("TCGATCGA", "8M", "8", 1, "chr1", Some(Seq(10, 15, 20, 25, 10, 15, 20, 25))),
      TestUtil.makeRead("TCGATCGA", "8M", "8", 1, "chr1", Some(Seq(10, 15, 20, 25, 10, 15, 20, 25))),
      TestUtil.makeRead("TCGACCCTCGA", "4M3I4M", "8", 1, "chr1", Some(Seq(10, 15, 20, 25, 5, 5, 5, 10, 15, 20, 25))))

    val lastPileup = Pileup(reads, "chr1", 8)
    lastPileup.elements.foreach(e => assertBases(e.sequencedBases, "A"))
    lastPileup.elements.forall(_.sequencedBases.headOption.exists(_ == Bases.A)) should be(true)

    lastPileup.elements.forall(_.isMatch) should be(true)
    lastPileup.elements.forall(_.qualityScore == 25) should be(true)
  }

  sparkTest("Load pileup from SAM file") {
    val pileup = loadPileup("same_start_reads.sam", 0)
    pileup.elements.length should be(10)
  }

  sparkTest("First 60 loci should have all 10 reads") {
    val pileup = loadPileup("same_start_reads.sam", 0)
    for (i <- 1 to 59) {
      val nextPileup = pileup.atGreaterLocus(i, Bases.N, Seq.empty.iterator)
      nextPileup.elements.length should be(10)
    }
  }

  sparkTest("test pileup element creation") {
    val read = TestUtil.makeRead("AATTG", "5M", "5", 0, "chr1")
    val firstElement = pileupElementFromRead(read, 0)

    firstElement.isMatch should be(true)
    firstElement.indexWithinCigarElement should be(0L)

    val secondElement = advancePileupElement(firstElement, 1L)
    secondElement.isMatch should be(true)
    secondElement.indexWithinCigarElement should be(1L)

    val thirdElement = advancePileupElement(secondElement, 2L)
    thirdElement.isMatch should be(true)
    thirdElement.indexWithinCigarElement should be(2L)

  }

  sparkTest("test pileup element creation with multiple cigar elements") {
    val read = TestUtil.makeRead("AAATTT", "3M3M", "6", 0, "chr1")

    val secondMatch = pileupElementFromRead(read, 3)
    secondMatch.isMatch should be(true)
    secondMatch.indexWithinCigarElement should be(0L)

    val secondMatchSecondElement = pileupElementFromRead(read, 4)
    secondMatchSecondElement.isMatch should be(true)
    secondMatchSecondElement.indexWithinCigarElement should be(1L)

  }

  test("insertion at contig start includes trailing base") {
    val contigStartInsertionRead = TestUtil.makeRead("AAAAAACGT", "5I4M", "4", 0, "chr1")
    val pileup = pileupElementFromRead(contigStartInsertionRead, 0)
    pileup.alignment should equal(Insertion("AAAAAA", List(31, 31, 31, 31, 31, 31)))
  }

  test("pileup alignment at insertion cigar-element throws") {
    val contigStartInsertionRead = TestUtil.makeRead("AAAAAACGT", "5I4M", "4", 0, "chr1")
    val pileup = PileupElement(
      read = contigStartInsertionRead,
      locus = 1,
      referenceBase = Bases.N,
      readPosition = 0,
      cigarElementIndex = 0,
      cigarElementLocus = 1,
      indexWithinCigarElement = 0
    )
    the[InvalidCigarElementException] thrownBy pileup.alignment
  }

  sparkTest("test pileup element creation with deletion cigar elements") {
    val read = TestUtil.makeRead("AATTGAATTG", "5M1D5M", "5^C5", 0, "chr1")
    val firstElement = pileupElementFromRead(read, 0)

    firstElement.isMatch should be(true)
    firstElement.indexWithinCigarElement should be(0L)

    val deletionElement = advancePileupElement(firstElement, 4L)
    deletionElement.alignment should equal(Deletion("GC", 1.toByte))
    deletionElement.isDeletion should be(true)
    deletionElement.indexWithinCigarElement should be(4L)

    val midDeletionElement = advancePileupElement(deletionElement, 5L)
    midDeletionElement.isMidDeletion should be(true)
    midDeletionElement.indexWithinCigarElement should be(0L)

    val pastDeletionElement = advancePileupElement(midDeletionElement, 6L)
    pastDeletionElement.isMatch should be(true)
    pastDeletionElement.indexWithinCigarElement should be(0L)

    val continuePastDeletionElement = advancePileupElement(pastDeletionElement, 9L)
    continuePastDeletionElement.isMatch should be(true)
    continuePastDeletionElement.indexWithinCigarElement should be(3L)

  }

  sparkTest("Loci 10-19 deleted from half of the reads") {
    val pileup = loadPileup("same_start_reads.sam", 0)
    val deletionPileup = pileup.atGreaterLocus(9, Bases.A, Seq.empty.iterator)
    deletionPileup.elements.map(_.alignment).count {
      case Deletion(bases, _) => {
        Bases.basesToString(bases) should equal("AAAAAAAAAAA")
        true
      }
      case _ => false
    } should be(5)
    for (i <- 10 to 19) {
      val nextPileup = pileup.atGreaterLocus(i, Bases.N, Seq.empty.iterator)
      nextPileup.elements.count(_.isMidDeletion) should be(5)
    }
  }

  sparkTest("Loci 60-69 have 5 reads") {
    val pileup = loadPileup("same_start_reads.sam", 0)
    for (i <- 60 to 69) {
      val nextPileup = pileup.atGreaterLocus(i, Bases.N, Seq.empty.iterator)
      nextPileup.elements.length should be(5)
    }
  }

  sparkTest("Pileup.Element basic test") {
    intercept[NullPointerException] {
      val e = pileupElementFromRead(null, 42)
    }

    val decadentRead1 = testAdamRecords(0)

    // read1 starts at SAM:6 → 0-based 5
    // and has CIGAR: 29M10D31M
    // so, the length is 70
    intercept[AssertionError] {
      pileupElementFromRead(decadentRead1, 0)
    }
    intercept[AssertionError] {
      pileupElementFromRead(decadentRead1, 4)
    }
    intercept[AssertionError] {
      pileupElementFromRead(decadentRead1, 5 + 70)
    }
    val at5 = pileupElementFromRead(decadentRead1, 5)
    assert(at5 != null)
    assertBases(at5.sequencedBases, "A")
    assert(at5.sequencedBases.headOption.exists(_ == Bases.A))

    // At the end of the read:
    assert(pileupElementFromRead(decadentRead1, 74) != null)
    intercept[AssertionError] {
      pileupElementFromRead(decadentRead1, 75)
    }

    // Just before the deletion
    val deletionPileup = pileupElementFromRead(decadentRead1, 5 + 28)
    deletionPileup.alignment should equal(Deletion("AGGGGGGGGGG", 1.toByte))

    // Inside the deletion
    val at29 = pileupElementFromRead(decadentRead1, 5 + 29)
    assert(at29.sequencedBases.size === 0)
    val at38 = pileupElementFromRead(decadentRead1, 5 + 38)
    assert(at38.sequencedBases.size === 0)
    // Just after the deletion
    assertBases(pileupElementFromRead(decadentRead1, 5 + 39).sequencedBases, "A")

    // advanceToLocus is a no-op on the same locus,
    // and fails in lower loci
    forAll(Table("locus", List(5, 33, 34, 43, 44, 74): _*)) { locus =>
      val elt = pileupElementFromRead(decadentRead1, locus)
      assert(advancePileupElement(elt, locus) === elt)
      intercept[AssertionError] {
        advancePileupElement(elt, locus - 1)
      }
      intercept[AssertionError] {
        advancePileupElement(elt, 75)
      }
    }

    val read3Record = testAdamRecords(2) // read3
    val read3At15 = pileupElementFromRead(read3Record, 15)
    assert(read3At15 != null)
    assertBases(read3At15.sequencedBases, "A")
    assertBases(advancePileupElement(read3At15, 16).sequencedBases, "T")
    assertBases(advancePileupElement(read3At15, 17).sequencedBases, "C")
    assertBases(advancePileupElement(advancePileupElement(read3At15, 16), 17).sequencedBases, "C")
    assertBases(advancePileupElement(read3At15, 18).sequencedBases, "G")
  }

  sparkTest("Read4 has CIGAR: 10M10I10D40M; ACGT repeated 15 times") {
    // Read4 has CIGAR: 10M10I10D40M
    // It's ACGT repeated 15 times
    val decadentRead4 = testAdamRecords(3)
    val read4At20 = pileupElementFromRead(decadentRead4, 20)
    assert(read4At20 != null)
    for (i <- 0 until 2) {
      assert(advancePileupElement(read4At20, 20 + i * 4 + 0).sequencedBases(0) == 'A')
      assert(advancePileupElement(read4At20, 20 + i * 4 + 1).sequencedBases(0) == 'C')
      assert(advancePileupElement(read4At20, 20 + i * 4 + 2).sequencedBases(0) == 'G')
      assert(advancePileupElement(read4At20, 20 + i * 4 + 3).sequencedBases(0) == 'T')
    }

    val read4At30 = advancePileupElement(read4At20, 20 + 9)
    read4At30.isInsertion should be(true)
    (read4At30.sequencedBases: String) should equal("CGTACGTACGT")
  }

  sparkTest("Read5: ACGTACGTACGTACG, 5M4=1X5=") {

    // Read5: ACGTACGTACGTACG, 5M4=1X5=, [10; 25[
    //        MMMMM====G=====
    val decadentRead5 = testAdamRecords(4)
    val read5At10 = pileupElementFromRead(decadentRead5, 10)
    assert(read5At10 != null)
    assertBases(advancePileupElement(read5At10, 10).sequencedBases, "A")
    assertBases(advancePileupElement(read5At10, 14).sequencedBases, "A")
    assertBases(advancePileupElement(read5At10, 18).sequencedBases, "A")
    assertBases(advancePileupElement(read5At10, 19).sequencedBases, "C")
    assertBases(advancePileupElement(read5At10, 20).sequencedBases, "G")
    assertBases(advancePileupElement(read5At10, 21).sequencedBases, "T")
    assertBases(advancePileupElement(read5At10, 22).sequencedBases, "A")
    assertBases(advancePileupElement(read5At10, 24).sequencedBases, "G")
  }

  sparkTest("read6: ACGTACGTACGT 4=1N4=4S") {
    // Read6: ACGTACGTACGT 4=1N4=4S
    // one `N` and soft-clipping at the end
    val decadentRead6 = testAdamRecords(5)
    val read6At40 = pileupElementFromRead(decadentRead6, 40)
    assert(read6At40 != null)
    assertBases(advancePileupElement(read6At40, 40).sequencedBases, "A")
    assertBases(advancePileupElement(read6At40, 41).sequencedBases, "C")
    assertBases(advancePileupElement(read6At40, 42).sequencedBases, "G")
    assertBases(advancePileupElement(read6At40, 43).sequencedBases, "T")
    assertBases(advancePileupElement(read6At40, 44).sequencedBases, "")
    assertBases(advancePileupElement(read6At40, 45).sequencedBases, "A")
    assertBases(advancePileupElement(read6At40, 48).sequencedBases, "T")
    intercept[AssertionError] {
      advancePileupElement(read6At40, 49).sequencedBases
    }
  }

  sparkTest("read7: ACGTACGT 4=1N4=4H, one `N` and hard-clipping at the end") {
    val decadentRead7 = testAdamRecords(6)
    val read7At40 = pileupElementFromRead(decadentRead7, 40)
    assert(read7At40 != null)
    assertBases(advancePileupElement(read7At40, 40).sequencedBases, "A")
    assertBases(advancePileupElement(read7At40, 41).sequencedBases, "C")
    assertBases(advancePileupElement(read7At40, 42).sequencedBases, "G")
    assertBases(advancePileupElement(read7At40, 43).sequencedBases, "T")
    assertBases(advancePileupElement(read7At40, 44).sequencedBases, "")
    assertBases(advancePileupElement(read7At40, 45).sequencedBases, "A")
    assertBases(advancePileupElement(read7At40, 48).sequencedBases, "T")
    intercept[AssertionError] {
      advancePileupElement(read7At40, 49).sequencedBases
    }
  }

  test("create and advance pileup element from RNA read") {
    val rnaRead = TestUtil.makeRead(
      sequence = "CCCCAGCCTAGGCCTTCGACACTGGGGGGCTGAGGGAAGGGGCACCTGCC",
      cigar = "7M191084N43M",
      mdtag = "9T24T7G7",
      start = 229538779,
      chr = "chr1")

    val rnaPileupElement = PileupElement(rnaRead, 229538779L, Bases.C)

    // Second base
    assertBases(advancePileupElement(rnaPileupElement, 229538780L).sequencedBases, "C")

    // Third base
    assertBases(advancePileupElement(rnaPileupElement, 229538781L).sequencedBases, "C")

    // In intron
    assertBases(advancePileupElement(rnaPileupElement, 229539779L).sequencedBases, "")

    // Last base
    assertBases(advancePileupElement(rnaPileupElement, 229729912L).sequencedBases, "C")

  }

  sparkTest("create pileup from RNA reads") {
    val rnaReadsPileup = loadPileup("testrna.sam", locus = 229580594)

    // 94 reads in the testrna.sam
    // 3 reads end at 229580707 and 1 extends further
    rnaReadsPileup.depth should be(94)

    val movedRnaReadsPileup = rnaReadsPileup.atGreaterLocus(229580706, Bases.A, Iterator.empty)
    movedRnaReadsPileup.depth should be(4)

    movedRnaReadsPileup.atGreaterLocus(229580707, Bases.N, Iterator.empty).depth should be(1)
  }

  test("pileup in the middle of a deletion") {
    val reads = Seq(
      TestUtil.makeRead("TCGAAAAGCT", "5M6D5M", "5^GCTTCG5", 0),
      TestUtil.makeRead("TCGAAAAGCT", "5M6D5M", "5^GCTTCG5", 0),
      TestUtil.makeRead("TCGAAAAGCT", "5M6D5M", "5^GCTTCG5", 0)
    )
    val deletionPileup = Pileup(reads, "chr1", 4)
    val deletionAlleles = deletionPileup.distinctAlleles
    deletionAlleles.size should be(1)
    deletionAlleles(0) should be(Allele("AGCTTCG", "A"))

    val midDeletionPileup = Pileup(reads, "chr1", 5)
    val midAlleles = midDeletionPileup.distinctAlleles
    midAlleles.size should be(1)
    midAlleles(0) should be(Allele("G", ""))
  }

}

