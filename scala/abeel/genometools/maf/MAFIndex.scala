package abeel.genometools.maf

import abeel.genometools.Main
import java.io.File
import net.sf.jannot.source.Locator
import net.sf.jannot.mafix.MafixFactory
import java.io.FileInputStream
import java.io.FileInputStream

object MAFIndex extends Main {
  case class Config(val inputFile: File = null, val outputFile: File = null)

  override val description = "Tool to index maf files."

  override val version = """
    2016/10/26       Initial version included in genometools
   """

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar mafix") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input MAF file. ")
      opt[File]('o', "output") action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written, default is <input>.bgz and <input>.bgz.mfi")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)

      processFile(config)

    }
  }

  def processFile(config: Config) {
    val input = config.inputFile
    val output = if (config.outputFile != null) config.outputFile else new File(config.inputFile + ".bgz")
    val idx = new File(output + ".mfi")
    MafixFactory.generateBlockZippedFile(new FileInputStream(input), output);
    MafixFactory.generateIndex(new Locator(output).stream(), idx);
  }

}