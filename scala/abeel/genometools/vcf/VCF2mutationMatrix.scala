package abeel.genometools.vcf

import java.io.File
import java.io.PrintWriter
import atk.util.Tool
import abeel.genometools.Main
import abeel.genometools.GenomeToolsConsole
import atk.io.ExtensionFileFilter
import atk.compbio.vcf.VCFFile
import scala.collection.JavaConversions._
import atk.util.NaturalOrderComparator

object VCF2mutationMatrix extends Main {

  override val description = """ Tool to convert a set of VCF files to a mutation matrix

WARNING: THIS TOOL HAS HARDCODED MAGIC VALUES AND IS INCOMPLETE!!!

This tool is still in development and is not for general use.




"""

  case class Config(val input: File = null, val output: File = null, val sizeLimit: Int = 10)

  override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar vcf2matrix") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Input file with a list of VCF files")
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Output file")
      opt[Int]("sizeLimit") action { (x, c) => c.copy(sizeLimit = x) } text ("Size limit of mutations to consider. Sites in the genome that have variants larger than this limit will be output in a separate file. (Default = 10)")
      //      opt[Int]("window") action { (x, c) => c.copy(window= x) } text ("Window size")
      //      opt[Unit]("incomplete")action { (x, c) => c.copy(incomplete= true) } text ("Include incomplete windows")

    }
    parser.parse(args, Config()) map { config =>

      init(config.input + ".log")
      val pw = new PrintWriter(config.output)
      pw.println(generatorInfo(config))
      pw.println()

      var time = System.currentTimeMillis()

      val files =tLines(config.input).map(new File(_))
      pw.println("# Input files")
      pw.println("# " + files.mkString("\n# "))

      val pxB = new PrintWriter(config.output + ".ignored")
      pxB.println(generatorInfo(config))

      val blackListPositions = for (file <- files) yield {
        file.getName -> VCFFile(file).filter(vcf => vcf.refLength > config.sizeLimit).map(vcf => vcf.pos -> (vcf.pos + vcf.refLength)).toList

      }
      pxB.println(blackListPositions.map(f => f._1 + "\t" + f._2.mkString(", ")).mkString("\n"))
      pxB.close()

      val vcfs = (for (file <- files) yield {
        file.getName() -> VCFFile(file).filter(vcf => vcf.refLength <= config.sizeLimit && vcf.altLength <= config.sizeLimit).map(vcf => vcf.pos + "_" + vcf.ref + "_" + vcf.alt -> vcf.pass).partition(_._2)

      }).toMap

      val pass = vcfs.mapValues(_._1.map(_._1).toList)
      val amb = vcfs.mapValues(_._2.map(_._1).toList)

      val variantSet = pass.values.flatten.toSet.toList.sorted(naturalOrdering)

      println(variantSet.mkString("\n"))

      //      val file = VCFFile(config.input)
      //      val linesIn = file.toList
      //      println("read: " + new TimeInterval(System.currentTimeMillis() - time))
      //      time = System.currentTimeMillis()

      pw.close()
      finish()

    }
  }
}
