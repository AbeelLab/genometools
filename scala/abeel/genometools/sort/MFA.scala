package abeel.genometools.sort

import java.io.File
import abeel.genometools.Main
import atk.util.NaturalOrderComparator
import net.sf.jannot.source.FileSource
import net.sf.jannot.parser.FastaParser
import scala.collection.JavaConversions._
import java.io.PrintWriter

object MFA extends Main {

  override val description = """Tool to sort mfa files"""

  case class Config(val inputFile: File = null, val outputFile: File = null, val n:Int=80)

  override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar sort_mfa") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input mfa formatted file.")
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("Output mfa formatted file with entries sorted in alphanumerci order.")
      opt[Int]('n', "n") action { (x, c) => c.copy(n = x) } text ("Characters per line in mfa")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)

      FastaParser.forceEntries=true
      
      val es = new FileSource(config.inputFile).read()

      val map = es.map{e =>
        
        (e.getID() -> e.sequence().stringRepresentation())}.toMap

      //      val mfa=tLines(config.inputFile)

      //      val map=mfa.grouped(2).map(f=>f(0)->f(1)).toMap

      val pw = new PrintWriter(config.outputFile)

      val sorted = map.keySet.toList.sorted(naturalOrdering)

     
      sorted.map { k =>
        pw.println(">" + k)
        pw.println(map(k).grouped(config.n).mkString("\n"))

      }

      pw.close

    }

  }

}