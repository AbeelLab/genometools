package abeel.genometools.faq

import java.util.Properties
import java.io.File
import atk.util.Tool
import net.sf.samtools.SAMFileReader
import scala.collection.JavaConversions._
import be.abeel.util.FrequencyMap
import be.abeel.util.FrequencyMapUtils
import abeel.genometools.Main
import java.util.HashMap
import atk.compbio.fastq.FastQFile
import atk.compbio.DNAString
import atk.compbio.fastq.FastAFile
import java.io.PrintWriter

object FaqStats extends Main {

  case class Config(val inputFile: File = null, val outputFile: File = null, val fastq: Boolean = false)

  override val description = "Tool to calculate statistics from fasta or fastq file."

  override val version = """
    2016/09/12       Initial version included in genometools
   """

  override def main(args: Array[String]) {

   
    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar faqstats") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input file. By default FASTA formatted. If you have a FASTQ, use the --fq flag") //, { v: String => config.spacerFile = v })
      opt[Unit]("fq") action { (x, c) => c.copy(fastq = true) } text ("If you have a FASTQ file, use this flag") //, { v: String => config.spacerFile = v })
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)

      processFile(config)

    }
  }

  private def processFile(config: Config) {

    val seqIt =
      if (config.fastq) {
        FastQFile(config.inputFile).map(_.seq)
      } else {
        FastAFile(config.inputFile).map(_.seq)
      }

    val ds = seqIt.map(new DNAString(_))

    val map = scala.collection.mutable.Map[String, Int]().withDefaultValue(0)
    for (seq <- ds) {

      if (seq.hasN)
        map("n") += 1
      if (seq.hasGap)
        map("gap") += 1
      if (seq.hasUndefined)
        map("undefined") += 1

      map("readLen_" + seq.size) += 1
      map("count") += 1
    }

    val pw = new PrintWriter(config.outputFile)
    pw.println(generatorInfo(config))
    map.map { case (x, y) => pw.println(x + "\t" + y) }

    pw.close

    pw.close

  }

}