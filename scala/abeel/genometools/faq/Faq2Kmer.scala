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

object Faq2Kmer extends Main {

  case class Config(val inputFile: File = null, val outputFile: File = null, val fastq: Boolean = false, val kmer: Int = 4)

  override val description = "Tool to draw kmer statistics from fasta or fastq file."

  override val version = """
    2016/09/12       Initial version included in genometools
   """

  override def main(args: Array[String]) {

   

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar faq2kmer") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input file. By default FASTA formatted. If you have a FASTQ, use the --fq flag")
      opt[Unit]("fq") action { (x, c) => c.copy(fastq = true) } text ("If you have a FASTQ file, use this flag") 
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")
      opt[Int]('k', "kmer") action { (x, c) => c.copy(kmer = x) } text ("Kmer length, default = 4")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)

      processFile(config)

    }
  }

  private def processFile(config: Config) {

    val kmerIt = if (config.fastq)
      FastQFile(config.inputFile).map(fr => fr.seq).map(seq => seq.sliding(config.kmer))
    else
      FastAFile(config.inputFile).map(fr => fr.seq).map(seq => seq.sliding(config.kmer))

    val map = scala.collection.mutable.Map[DNAString, Int]().withDefaultValue(0)
    for (sq <- kmerIt.flatten) {
      val seq = new DNAString(sq)
      map(seq) += 1

    }

    val pw = new PrintWriter(config.outputFile)
    pw.println(generatorInfo(config))
    map.map { case (x, y) => pw.println(x + "\t" + y) }

    pw.close

  }

}