package abeel.genometools.vcf

import atk.util.Tool
import java.io.File
import java.io.PrintWriter
import atk.compbio.vcf.VCFFile
import atk.compbio.vcf.VCFLine
import atk.util.MD5Tools
import scala.io.Source

object VCFDiff extends Tool {

  def main(args: Array[String]): Unit = {

    case class Config(val aOut: File = null, val bOut: File = null, val a: File = null, val b: File = null)

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar vcf-diff") {
      opt[File]("output-a") action { (x, c) => c.copy(aOut = x) } text ("File with unique mutations in first input file")
      opt[File]("output-b") action { (x, c) => c.copy(bOut = x) } text ("File with unique mutations in second input file")
      opt[File]('a', "first") required () action { (x, c) => c.copy(a = x) } text ("First VCF file you want to analyze")
      opt[File]('b', "second") required () action { (x, c) => c.copy(b = x) } text ("Second VCF file you want to analyze")
    }

    parser.parse(args, Config()) map { config =>

      /* Build diff set */
      val a = VCFFile(config.a)
      val b = VCFFile(config.b)

      val f: VCFLine => Boolean = {
        f =>
          !f.ref.equals(f.alt)
      }

      val df: VCFLine => String = { f =>

        f.refGenome + "." + f.pos + "." + f.ref + "." + f.alt

      }

      val fa = a.filter(f).map(df)
      val fb = a.filter(f).map(df)

      val diff = fa.toSet.diff(fb.toSet)

      def process(xFile: File, diff: Set[String], xOutput: File) = {

        val pw = new PrintWriter(xOutput)
        val header=Source.fromFile(xFile).getLines.filter(_.startsWith("#"))
        
        val x = VCFFile(xFile)
        for (line <- x) {
          val id = df(line)
          if (diff.contains(id))
            pw.println(line.line)
        }
        pw.println("## This file is complete!")
        pw.close

      }

      process(config.a, diff, config.aOut)
      process(config.a, diff, config.bOut)

    }

  }

}
