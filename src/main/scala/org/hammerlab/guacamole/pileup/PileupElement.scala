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

import htsjdk.samtools.{ CigarElement, CigarOperator }
import org.hammerlab.guacamole.{ Bases, CigarUtils }
import org.hammerlab.guacamole.reads.MappedRead
import org.hammerlab.guacamole.variants.Allele

import scala.annotation.tailrec

/**
 * A [[PileupElement]] represents the bases sequenced by a particular read at a particular reference locus.
 *
 * @param read The read this [[PileupElement]] is coming from.
 * @param locus The reference locus.
 * @param readPosition The offset into the sequence of bases in the read that this element corresponds to.
 * @param cigarElementIndex The index in the read's sequence of cigar elements ([[org.hammerlab.guacamole.reads.MappedRead.cigarElements]])
 *                          of the element that contains the current readPosition.
 * @param cigarElementLocus The reference START position of the current cigar element.
 *                          If the element is an INSERTION this the PRECEDING reference base
 * @param indexWithinCigarElement The offset of this element within the current cigar element.
 */
case class PileupElement(
    read: MappedRead,
    locus: Long,
    referenceBase: Byte,
    readPosition: Int,
    cigarElementIndex: Int,
    cigarElementLocus: Long,
    indexWithinCigarElement: Int) {

  assume(locus >= read.start)
  assume(locus < read.end)

  lazy val cigarElement = read.cigarElements(cigarElementIndex)
  lazy val nextCigarElement =
    if (cigarElementIndex + 1 < read.cigarElements.size) {
      Some(read.cigarElements(cigarElementIndex + 1))
    } else {
      None
    }

  lazy val referenceStringIdx =
    (cigarElementLocus - read.start).toInt +
      (if (cigarElement.getOperator.consumesReferenceBases()) indexWithinCigarElement else 0)

  lazy val cigarElementReadLength = CigarUtils.getReadLength(cigarElement)
  lazy val cigarElementReferenceLength = CigarUtils.getReferenceLength(cigarElement)
  lazy val cigarElementEndLocus = cigarElementLocus + cigarElementReferenceLength

  lazy val alignment: Alignment = {

    /*
     * True if this is the last base of the current cigar element.
     */
    def isFinalCigarBase: Boolean = indexWithinCigarElement == cigarElement.getLength - 1
    val cigarOperator = cigarElement.getOperator
    val nextBaseCigarElement = if (isFinalCigarBase) nextCigarElement else Some(cigarElement)
    val nextBaseCigarOperator = nextBaseCigarElement.map(_.getOperator)

    def makeInsertion(cigarElem: CigarElement) =
      Insertion(
        read.sequence.view(
          readPosition,
          readPosition + CigarUtils.getReadLength(cigarElem) + 1
        ),
        read.baseQualities.view(
          readPosition,
          readPosition + CigarUtils.getReadLength(cigarElem) + 1
        )
      )

    (cigarOperator, nextBaseCigarOperator) match {

      // Since insertions by definition have no corresponding reference loci, there is a choice in whether we "attach"
      // them to the preceding or following locus. Here we attach them to the preceding base, since that seems to be the
      // conventional choice. That is, if we have a match followed by an insertion, the final base of the match will
      // get combined with the insertion into one Alignment, at the match's reference locus.
      case (CigarOperator.M, Some(CigarOperator.I)) | (CigarOperator.EQ, Some(CigarOperator.I)) =>
        makeInsertion(nextCigarElement.get)

      // The exception to the above is insertion at the start of a contig, where there is no preceding reference base to
      // anchor to; in this case, the spec calls for including the reference base immediately after the insertion (the
      // first reference base of the contig).
      case (CigarOperator.I, Some(_)) if cigarElementLocus == 0 =>
        makeInsertion(cigarElement)

      // In general, a PileupElement pointing at an Insertion cigar-element is an error.
      case (CigarOperator.I, _) => throw new InvalidCigarElementException(this)

      case (CigarOperator.M | CigarOperator.EQ | CigarOperator.X, Some(CigarOperator.D)) =>
        val deletedBases =
          referenceBase ::
            (referenceStringIdx + 1 until referenceStringIdx + 1 + nextCigarElement.get.getLength)
            .map(offset => read.mdTagOpt.get.deletions(read.start + offset).toByte).toList
        val anchorBaseSequenceQuality = read.baseQualities(readPosition)
        Deletion(deletedBases, anchorBaseSequenceQuality)
      case (CigarOperator.D, _) =>
        // TODO(ryan): can a cigar begin with a 'D' operator?
        MidDeletion(read.mdTagOpt.get.deletions(locus).toByte)
      case (op, Some(CigarOperator.D)) =>
        // TODO(ryan): are there sane cases where a 'D' is not preceded by an 'M'?
        throw new AssertionError(
          "Found deletion preceded by cigar operator %s at PileupElement for read %s at locus %d".format(op, read.toString, locus)
        )
      case (CigarOperator.M, _) | (CigarOperator.EQ, _) | (CigarOperator.X, _) =>
        val base: Byte = read.sequence(readPosition)
        val quality = read.baseQualities(readPosition)
        if (base == referenceBase) {
          Match(base, quality)
        } else {
          Mismatch(base, quality, referenceBase)
        }
      case (CigarOperator.S, _) | (CigarOperator.N, _) | (CigarOperator.H, _) => Clipped
      case (CigarOperator.P, _) =>
        throw new AssertionError("`P` CIGAR-ops should have been ignored earlier in `findNextCigarElement`")
    }
  }

  /* If you only care about what kind of CigarOperator is at this position, but not its associated sequence, then you
   * can use these state variables.
   */
  lazy val isInsertion = alignment match { case Insertion(_, _) => true; case _ => false }
  lazy val isDeletion = alignment match { case Deletion(_, _) => true; case _ => false }
  lazy val isMidDeletion = alignment match { case MidDeletion(_) => true; case _ => false }
  lazy val isMismatch = alignment match { case Mismatch(_, _, _) => true; case _ => false }
  lazy val isMatch = alignment match { case Match(_, _) => true; case _ => false }

  /**
   * The sequenced nucleotides at this element.
   *
   * If the current element is a deletion, then this is the empty array. If it's
   * an insertion, then this will be an array of length >= 1: the contents of
   * the inserted sequence starting at the current locus. Otherwise, this is
   * an array of length 1.
   */
  lazy val sequencedBases: Seq[Byte] = alignment.sequencedBases
  lazy val referenceBases: Seq[Byte] = alignment.referenceBases

  lazy val allele: Allele = Allele(referenceBases, sequencedBases)

  /*
   * Base quality score, phred-scaled.
   *
   * For matches and mismatches this is the base quality score of the current base.
   * For insertions this the minimum base quality score of all the bases in the insertion.
   * For deletions this is the mapping quality as there are no base quality scores available.
   */
  lazy val qualityScore: Int = alignment match {
    case Clipped | MidDeletion(_) => read.alignmentQuality
    case Deletion(_, qs)          => qs
    case MatchOrMisMatch(_, qs)   => qs
    case Insertion(_, qss)        => qss.min
  }

  /**
   * Returns a new [[PileupElement]] of the same read, advanced by one cigar element.
   */
  def advanceToNextCigarElement: PileupElement = {
    val readPositionOffset =
      if (cigarElement.getOperator.consumesReadBases()) {
        // If this [[CigarElement]] consumes read bases then [[readPosition]] will advance by the rest of this
        // [[CigarElement]].
        cigarElement.getLength - indexWithinCigarElement
      } else {
        // Otherwise, [[readPosition]] will stay the same.
        0
      }

    val nextLocus = locus + (cigarElementReferenceLength - indexWithinCigarElement)
    PileupElement(
      read,
      nextLocus,
      Bases.N, // placeholder before we advance to a proper locus
      readPosition + readPositionOffset,
      cigarElementIndex + 1,
      cigarElementLocus + cigarElementReferenceLength,
      // Even if we are somewhere in the middle of the current cigar element, lock to the beginning of the next one.
      indexWithinCigarElement = 0
    )
  }

  /**
   * Returns whether the current cigar element of this [[org.hammerlab.guacamole.reads.MappedRead]] contains the given reference locus.
   *
   * Can only return true if the cigar element consumes reference bases.
   */
  def currentCigarElementContainsLocus(referenceLocus: Long): Boolean = {
    cigarElementLocus <= referenceLocus && referenceLocus < cigarElementEndLocus
  }

  /**
   * Returns a new [[PileupElement]] of the same read at a different locus.
   *
   * To enable an efficient implementation, newLocus must be greater than the current locus.
   *
   * @param newLocus The desired locus of the new [[PileupElement]]. It must be greater than the current locus, and
   *                 not past the end of the current read.
   * @param newReferenceBase The reference base at [[newLocus]]
   *
   * @return A new [[PileupElement]] at the given locus.
   */
  @tailrec
  final def advanceToLocus(newLocus: Long, newReferenceBase: Byte): PileupElement = {
    assume(newLocus >= locus, "Can't rewind to locus %d from %d. Pileups only advance.".format(newLocus, locus))
    assume(newLocus < read.end, "This read stops at position %d. Can't advance to %d".format(read.end, newLocus))
    if (currentCigarElementContainsLocus(newLocus)) {
      // Aside: the current cigar element must consume reference bases if we've gotten here.
      val readPositionOffset =
        if (cigarElement.getOperator.consumesReadBases()) {
          (newLocus - cigarElementLocus - indexWithinCigarElement).toInt
        } else {
          // If this cigar doesn't consume read bases, then advancing within it will not consume [[readPosition]]
          0
        }

      this.copy(
        locus = newLocus,
        referenceBase = newReferenceBase,
        readPosition = readPosition + readPositionOffset,
        indexWithinCigarElement = (newLocus - cigarElementLocus).toInt
      )
    } else if (newLocus == 0 && cigarElement.getOperator == CigarOperator.I) {
      // NOTE(ryan): this is the rare case where we allow a [[PileupElement]] to exist at a non-reference-consuming
      // CigarElement (namely, an Insertion at the start of a contig). It is correct for us to emit an Insertion
      // alignment in such a case, where typically an Insertion [[CigarElement]] would get skipped over by subsequent
      // [[advanceToLocus]] calls, since it represents a zero-width interval of reference bases.
      this
    } else
      advanceToNextCigarElement.advanceToLocus(newLocus, newReferenceBase)
  }

  /**
   * Distance from the end of the reading frame
   * If the read was positive then the sequencing end also corresponds to the read end positions
   * If the read was negative, the sequencing end is the mapped start position
   */
  lazy val distanceFromSequencingEnd = if (read.isPositiveStrand) read.end - locus else locus - read.start

}

object PileupElement {

  /**
   * Create a new [[PileupElement]] backed by the given read at the specified locus. The read must overlap the locus.
   */
  def apply(read: MappedRead, locus: Long, referenceBase: Byte): PileupElement = {
    PileupElement(
      read = read,
      locus = read.start,
      referenceBase = Bases.N,
      readPosition = 0,
      cigarElementIndex = 0,
      cigarElementLocus = read.start,
      indexWithinCigarElement = 0
    ).advanceToLocus(locus, referenceBase)
  }
}

case class InvalidCigarElementException(elem: PileupElement)
  extends Exception(
    "Should not have a PileupElement at non-reference-consuming cigar-operator I. " +
      "Locus: %d, readPosition: %d, cigar: %s (elem idx %d)".format(
        elem.locus,
        elem.readPosition,
        elem.read.cigar.toString, elem.cigarElementIndex
      )
  )
