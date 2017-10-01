package abeel.genometools.vcfstats

import java.io.File
import java.io.PrintWriter
import atk.util.Tool
import atk.compbio.vcf.VCFFile
import be.abeel.util.CountMap
import atk.compbio.vcf.VCFFile
import scala.collection.JavaConversions._
import abeel.genometools.Main
object VCFStatistics extends Tool with Main {
  case class Config(val input: File = null, val output: File = null)

  override def main(args: Array[String]): Unit = {
    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar vcf-statistics ") {
      opt[File]('o', "output") action { (x, c) => c.copy(output = x) } text ("File where you want the output to be written, by default the output is written to the console.")
      opt[File]('i', "input") required() action { (x, c) => c.copy(input = x) } text ("VCF file you want to analyze")
    }

    parser.parse(args, Config()) map { config =>
      val pw = if (config.output != null) new PrintWriter(config.output) else new PrintWriter(System.out)

      pw.println(generatorInfo)

      pw.println("input\t" + config.input)
      pw.println("output\t" + config.output)

      val file = VCFFile(config.input)

      val fmAll = new CountMap[String];
      val fmPass=new CountMap[String];
      
      for (line <- file) {
        fmAll.count(line.variation.strType)
        if(line.pass)
          fmPass.count(line.variation.strType)
      }
      pw.println("-- All events")
      pw.println("total\t"+fmAll.totalCount())
      pw.println(fmAll.map(f => f._1 + "\t" + f._2).mkString("\n"))
      pw.println("-- Passing filter events")
      pw.println("total\t"+fmPass.totalCount())
      pw.println(fmPass.map(f => f._1 + "\t" + f._2).mkString("\n"))
      pw.close
    }
  }

}