package abeel.genometools.kmer

import java.io.File
import abeel.genometools.Main
import atk.compbio.DNAHash
import atk.util.BitSetTools
import java.io.PrintWriter
import atk.util.TimeInterval
import java.time.LocalDateTime
import atk.io.NixWriter

object GroupAP extends Main {
  case class Config(val output: File = null, val input: File = null)

  override val description = "Group profiles from an AP matrix into a list of rows per profile."

  override val version = """
    2016/09/28       Initial version included in genometools
   """

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar groupprofile") {

      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("File where you want the output to be written")
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("File containing a AP matrix for each Kmer per strain")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.output != null)
      assume(config.input != null)

      processFile(config)

    }

  }

  private def processFile(config: Config) {

    val it = tLinesIterator(config.input)
    val files = it.next()

    var profileLen = 0
    var hashLen = 0

    val everything = for (line <- it) yield {
      
      progress(100000)

      val arr = line.split("\t")
      val hash = DNAHash.hash(arr(0))
      if (hashLen == 0)
        hashLen = arr(0).size
      val profile = BitSetTools.fromString(arr(1))
      if (profileLen == 0)
        profileLen = arr(1).size
      hash -> profile
    }
    
    
    val groups = everything.toList.groupBy(_._2).mapValues(_.map(_._1))
    
    println("Writing output...")
    val pw=new NixWriter(config.output,config)
    pw.println("# Unique profiles: " + groups.size)
    pw.println(files)
    for (group <- groups) {
      pw.println(BitSetTools.toString(group._1, profileLen)+"\t"+group._2.size + "\t" + (group._2.map { x => DNAHash.unhash(x, hashLen) }).mkString(";"))

    }
    println("Done!")
    pw.close

  }

}