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

package org.hammerlab.guacamole.filters

import org.apache.commons.math3.util.ArithmeticUtils

object FishersExactTest {

  /* Fischer's exact test, returned as a probability. */
  def apply(totalA: Int, totalB: Int, conditionA: Int, conditionB: Int): Double = {
    math.exp(ArithmeticUtils.binomialCoefficientLog(totalA, conditionA) +
      ArithmeticUtils.binomialCoefficientLog(totalB, conditionB) -
      ArithmeticUtils.binomialCoefficientLog(totalA + totalB, conditionA + conditionB))
  }

  /* Fischer's exact test, returned as log base 10 probability. */
  def asLog10(totalA: Int, totalB: Int, conditionA: Int, conditionB: Int): Double = {
    (ArithmeticUtils.binomialCoefficientLog(totalA, conditionA) +
      ArithmeticUtils.binomialCoefficientLog(totalB, conditionB) -
      ArithmeticUtils.binomialCoefficientLog(totalA + totalB, conditionA + conditionB)
      / Math.log(10))
  }

  def log10BinomialTestPValue(n: Int, k: Int, p: Double = 0.5): Double = {
    val q = 1.0 - p
    val logEResult = ArithmeticUtils.binomialCoefficientLog(n, k) + k * Math.log(p) + (n - k) * Math.log(q)
    logEResult / Math.log(10)
  }
}