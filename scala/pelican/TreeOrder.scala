package pelican

import atk.compbio.tree.Tree
import scala.collection.JavaConversions._
import java.io.File
import abeel.genometools.Main

object TreeOrder extends Main {

  case class Config(val tree: File = null)
  
  override def description = """List all nodes in the tree in the order they will appear in the figure."""
  
  override def main(args: Array[String]): Unit = {
    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar list-tree") {
      opt[File]('t', "tree") required () action { (x, c) => c.copy(tree = x) } text ("Phylogenetic tree file.")

    }

    parser.parse(args, Config()).map { config =>
      fix(config);
    }

  }

  def fix(config: Config): Unit = {

    val tree = new Tree(config.tree.toString())

    println(tree.getLeaves(tree.root).toList.map(_.getName()).mkString("\n"))

  }

}