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

package org.bdgenomics.guacamole

import com.esotericsoftware.kryo.Kryo
import org.bdgenomics.adam.serialization.{ ADAMKryoRegistrator }

class GuacamoleKryoRegistrator extends ADAMKryoRegistrator {
  override def registerClasses(kryo: Kryo) {
    // For debugging, it can be helpful to uncomment the following line. This will raise an error if we try to serialize
    // anything without a registered Kryo serializer.
    // kryo.setRegistrationRequired(true)

    // Register ADAM serializers.
    super.registerClasses(kryo)

    // Register Guacamole serializers.
    kryo.register(classOf[MappedRead], new MappedReadSerializer)
    kryo.register(classOf[UnmappedRead], new UnmappedReadSerializer)
    kryo.register(classOf[LociSet], new LociSetSerializer)
    kryo.register(classOf[LociSet.SingleContig], new LociSetSingleContigSerializer)
    kryo.register(classOf[LociMap[Long]], new LociMapLongSerializer)
    kryo.register(classOf[LociMap.SingleContig[Long]], new LociMapLongSingleContigSerializer)
  }
}