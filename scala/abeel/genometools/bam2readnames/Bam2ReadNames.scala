package abeel.genometools.bam2readnames

import java.util.Properties
import java.io.File
import atk.util.Tool
import java.io.PrintWriter
import net.sf.samtools.SAMFileReader
import scala.collection.JavaConversions._


object Bam2ReadNames extends Tool {

  case class Config(val inputFile: File = null, val outputFile: File = null)

  def main(args: Array[String]) {

    println("##----------------------------------------------")
    println("## Bam2ReadNames.scala")
    println("## ")
    println("## Tool to extract read names from a bam file")
    println("## ")
    println("## ")
    println("## The program will conclude with a message")
    println("## that indicates the run was successful.")
    println("## ")
    println("## By Thomas Abeel (tabeel@broadinstitute.org)")
    println("##----------------------------------------------")

    try {
      val prop = new Properties();
      prop.load(Bam2ReadNames.getClass().getResourceAsStream("/tool.properties"));
      println("## Program=" + prop.getProperty("program"));
      println("## Version=" + prop.getProperty("version"));
    } catch {
      case e: Exception =>
        System.err.println("Problem while reading version information");
        e.printStackTrace();
    }

    val parser = new scopt.OptionParser[Config]("java -jar bam2readnames.jar") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input file in which to inject data") //, { v: String => config.spacerFile = v })
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)

      processFile(config)

    }
  }

  private def processFile(config: Config) {

    

    val sam = new SAMFileReader(config.inputFile)

    val pw = new PrintWriter(config.outputFile)

    val names = sam.iterator().map(f => pw.println(f.getReadName()))
    while (names.hasNext) names.next

    pw.close
    sam.close()

  }

}