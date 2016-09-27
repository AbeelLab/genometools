package abeel.genometools.nwk

import java.io.File
import java.io.PrintWriter
import atk.compbio.tree.Tree
import scala.collection.JavaConversions._
import abeel.genometools.Main
import atk.compbio.tree.TreeNode

object Nwk2Nodes extends Main {

  override val description = """Tool to convert a nwk file to a list of nodes, each including the list of strains under that node."""

  case class Config(val input: File = null, val output: File = null)

  override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar nwk2nodes") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Input file")
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Output file")

    }
    parser.parse(args, Config()) map { config =>
      replace(config)

    }
  }

  def replace(config: Config) {

    val pw = new PrintWriter(config.output)
    pw.println(generatorInfo(config))
    val tree = new Tree(config.input.toString())

    traverse(tree.root, 0, 0)

    def traverse(node: TreeNode, lvl: Int, childIdx: Int): Unit = {
      if (!node.isLeaf()) {
        val leaves = tree.getLeaves(node).map(_.getName)
        pw.println("internal_" + lvl + "_" + childIdx +"\t"+leaves.size+ "\t" + leaves.mkString(";"))
        for (child <- node.children.zipWithIndex) {
          traverse(child._1, lvl + 1, child._2)
        }
      } else {
        pw.println("leaf_" + lvl + "_" + childIdx +"\t1\t" + node.getName)
      }

    }

    pw.close

  }
}