package org.hammerlab.guacamole.commands

import org.hammerlab.guacamole.commands.GermlineAssemblyCaller.Arguments
import org.hammerlab.guacamole.util.{GuacFunSuite, TestUtil}

class GermlineAssemblyCallerSuite extends GuacFunSuite {

  sparkTest("test assembly caller") {
    val input = "assemble-reads-set3-chr2-73613071.sam"


    val args = new Arguments
    args.reads = TestUtil.testDataPath(input)
    args.referenceFastaPath = TestUtil.testDataPath("chr2.fasta")
    args.parallelism = 2
    args.variantOutput = TestUtil.tmpFileName(suffix = ".vcf")
    args.snvWindowRange = 1000

    GermlineAssemblyCaller.Caller.run(args, sc)

  }

}
