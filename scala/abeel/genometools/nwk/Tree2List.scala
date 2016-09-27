package abeel.genometools.nwk

import java.io.File
import java.io.PrintWriter
import atk.compbio.tree.Tree
import scala.collection.JavaConversions._
import abeel.genometools.Main

object Tree2List extends Main{

  case class Config(val input: File = null, val output: File = null)

  override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar nwk2list") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Input file")
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Output file")

    }
    parser.parse(args, Config()) map { config =>
      replace(config)

    }
  }

  def replace(config: Config) {

    val pw=new PrintWriter(config.output)
    pw.println(generatorInfo(config))
    val tree=new Tree(config.input.toString())
    val leaves=tree.getLeaves(tree.root)
    leaves.map(l=>pw.println(l.getName))
    pw.close
    
  }
}