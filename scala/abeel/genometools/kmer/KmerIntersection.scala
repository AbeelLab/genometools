package abeel.genometools.kmer

import atk.util.Tool
import atk.compbio.DNAHash
import abeel.genometools.Main
import java.io.File
import java.io.PrintWriter
import atk.io.NixWriter
import java.time.LocalDateTime

object KmerIntersection extends Main {

  case class Config(val outputFile: File = null, val postfix: String = "", val strains: String = "")

  override val description = "Tool to reduce kmer file size by filtering by count."

  override val version = """
    2016/09/27       Initial version included in genometools
   """

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar intersect-kmer") {
      //      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input BAM file. ")
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")
      opt[String]("postfix") action { (x, c) => c.copy(postfix = x) } text ("String to append to identifiers to get kmer files.")
      opt[String]("strains") action { (x, c) => c.copy(strains = x) } text ("Semi-colon separated list of strain identifiers, e.g. strain1;strain2;strain3")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.outputFile != null)

      processFile(config)

    }


  }

  def processFile(config: Config) {

    val map = scala.collection.mutable.Map[Long, Int]().withDefaultValue(0)
    var hashSize = 0
    for (strain <- config.strains.split(";")) {
      println("processing: " + strain + "\t" + LocalDateTime.now())
      for (line <- tLinesIterator(strain + config.postfix)) {
        val seq = line.split("\t")(0)
        val ds = DNAHash.hash(seq) //new DNAString(seq)
        if (hashSize == 0)
          hashSize = seq.size
        map(ds) += 1
      }
    }
    val size = config.strains.split(";").size

    val pw = new PrintWriter(config.outputFile)
    pw.println(generatorInfo(config))

    val filtered = map.filter(_._2 == size)
    filtered.map(f => pw.println(DNAHash.unhash(f._1, hashSize) + "\t" + f._2))
    pw.println("# " + filtered.size + " shared kmers")
    pw.close
  }
}