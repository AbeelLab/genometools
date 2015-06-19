package abeel.genometools.phy2maf

import java.io.File
import abeel.genometools.Main
import java.io.PrintWriter

object Phy2Maf extends Main {

  override val description ="""Tool to convert phylip to multi-fasta file format."""
  
  case class Config(val inputFile: File = null, val outputFile: File = null, val n:Int=80)

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar phy2maf") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input phylip formatted file.")
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("Output file in maf format.")
      opt[Int]('n',"n") action { (x, c) => c.copy(n = x) } text ("Characters per line in mfa")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)

    
      val phy=tLines(config.inputFile)
      val pw=new PrintWriter(config.outputFile)
      phy.drop(1).map{line=>
        val arr=line.split("\\s+")
        assume(arr.length==2,"Unexpected number of elements on line: "+line)
        pw.println(">"+arr(0))
        pw.println(arr(1).grouped(config.n).mkString("\n"))
      }
      pw.close

    }
  }
}