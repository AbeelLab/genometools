package abeel.genometools.bam

import java.io.File
import net.sf.samtools.SAMFileReader
import atk.compbio.DNAString
import scala.collection.JavaConversions._
import java.io.PrintWriter
import abeel.genometools.Main
import net.sf.samtools.SAMFileReader.ValidationStringency
import atk.util.TimeInterval
import atk.compbio.DNAHash

object Bam2Kmer extends Main {
  case class Config(val inputFile: File = null, val outputFile: File = null, val kmer: Int = 4)

  override val description = "Tool to draw kmer statistics from a BAM file."

  override val version = """
    2016/09/26       Initial version included in genometools
   """

  override def main(args: Array[String]) {
 
    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar bam2kmer") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input BAM file. ")
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")
      opt[Int]('k', "kmer") action { (x, c) => c.copy(kmer = x) } text ("Kmer length, default = " + new Config().kmer)

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)

      processFile(config)

    }
  }

  private def processFile(config: Config) {

    val sfr = new SAMFileReader(config.inputFile)
    sfr.setValidationStringency(ValidationStringency.LENIENT)

    val map = scala.collection.mutable.Map[Long, Int]().withDefaultValue(0)
    var counter = 1
    var discard = 0
    val startTime = System.currentTimeMillis()
    for (samRecord <- sfr.iterator()) {
      if (counter % 10000 == 0) {
        val interval = System.currentTimeMillis() - startTime
        println("Processing: " + counter + "\t" + new TimeInterval(interval) + "\t" + nf.format(counter*1000L / (interval + .1)) + " reads/s")
      }
      counter += 1
      val sr = new String(samRecord.getReadBases)

      for (seq <- sr.sliding(config.kmer)) {
        try {
          val ds = DNAHash.hash(seq) //new DNAString(seq)
          map(ds) += 1
        } catch {
          case _: Throwable => discard +=1;
        }

      }

    }
    

    val pw = new PrintWriter(config.outputFile)
    pw.println(generatorInfo(config))
    pw.println("# Processed reads = "+counter)
    pw.println("# Discarded reads = "+discard)
    map.map { case (x, y) => pw.println(DNAHash.unhash(x,config.kmer) + "\t" + y) }

    pw.close

  }
}