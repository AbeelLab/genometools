package abeel.genometools.gfa
import abeel.genometools.Main
import java.io.File
import java.io.PrintWriter
import scala.collection.mutable.MutableList
import atk.tools.Histogram
import atk.tools.Histogram.HistogramConfig
import atk.io.NixWriter

/**
 * Program to calculate statistics on GFA files
 */

object GFAStatistics extends Main {

  override val description = """Tool to calculate basic statistics from a GFA file.

This tool is still in development and is not for general use.


"""

  case class Config(val inputFile: File = null, val outputFile: File = null)

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar gfa-statistics") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input GFA formatted file.")
      opt[File]('o', "output") action { (x, c) => c.copy(outputFile = x) } text ("Output file containing statistics.")

    }

    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)
      val pw = new NixWriter(config.outputFile, config)
      val cc = new Config(config.inputFile, null)

      val gfa = tLines(cc.inputFile)
      //      val pw = new PrintWriter(cc.outputFile)
      val map = gfa.map(line => line.split("\t"))

      println("Read in memory")
      val grouped = map.groupBy { x => x(0) }
      val headers = grouped("H")
      pw.println("#Headers: " + headers.flatten.toList)

      val genomeCount = headers(1)(1).split(";").size
      pw.println("Genome count: " + genomeCount)

      val segments = gfa.filter(_(0) == 'S').map(line => {
        val arr = line.split("\t")
        assume(arr(4).startsWith("ORI"))

        new Segment(arr(1).toInt, arr(2).length(), arr(4).drop(6).replaceAll(".fna", "").replaceAll(".fasta", "").split(";").toList.map(v => Dictionary.get(v)), MutableList.empty, MutableList.empty)
      })

      println("Mapping identifiers...")
      //map ID to segment
      val nodeMap = segments.map(segment => {
        segment.idx -> segment
      }).toMap

     

      /* parse all edges -- has side effect to update segments */
      val links = gfa.filter(_(0) == 'L').map(line => {
        val arr = line.split("\t")
        val from = arr(1).toInt
        val to = arr(3).toInt

        nodeMap(from).outgoing += to
        nodeMap(to).incoming += from

        from.toInt -> to.toInt

      })
      pw.println("")
      pw.println("Graph statistics")
      pw.println("----------------")
      val ss = grouped("S")
      val ll = grouped("L")
      pw.println("genomes=\t" + genomeCount)
      pw.println("segments=\t" + ss.size)
      pw.println("links=\t" + links.size)
      val sumLen = (ss.map(seg => seg(2).size)).sum
      pw.println("total segment size=\t" + sumLen)

      
      pw.println("")
      pw.println("Link statistics")
      pw.println("---------------")
      /* Incoming edge distribution */

      val inDegree = links.groupBy(_._2).mapValues(_.size)

      Histogram.plot(inDegree.toList.map(_._2 + 0.0), new HistogramConfig(nobin = true, x = "In-degree per segment", y = "Frequency", outputPrefix = config.outputFile + ".indegree", domainStart = 0, domainEnd = genomeCount))

      val histIn = inDegree.map(_._2).groupBy(identity).mapValues(_.size).toList.sortBy(_._1)
      pw.println()
      pw.println("In degree distribution: ")
      histIn.map { case (x, y) => pw.println(x + "\t" + y) }

      val outDegree = links.groupBy(_._1).mapValues(_.size)
      Histogram.plot(outDegree.toList.map(_._2 + 0.0), new HistogramConfig(nobin = true, x = "Out-degree per segment", y = "Frequency", outputPrefix = config.outputFile + ".outdegree", domainStart = 0, domainEnd = genomeCount))
      val histOut = outDegree.map(_._2).groupBy(identity).mapValues(_.size).toList.sortBy(_._1)
     
      pw.println()
      pw.println("Out degree distribution: ")
      histOut.map { case (x, y) => pw.println(x + "\t" + y) }

      pw.println("")
      pw.println("Node statistics")
      pw.println("---------------")

      /* Head nodes */
      val headNodes = segments.filter(seg => seg.incoming.size == 0)
      pw.println()
      pw.println("Root nodes: " + headNodes.size)
      pw.println(headNodes.mkString("\n"))

      val tailNodes = segments.filter(seg => seg.outgoing.size == 0)
      pw.println()
      pw.println("Tail nodes: " + tailNodes.size)
      pw.println(tailNodes.mkString("\n"))

      nfP.setMaximumFractionDigits(0)

      val genomeCountPerSegmentWithSegmentLength = (ss.map(seg => seg(4).split(";").size -> seg(2).size))

      pw.println()
      pw.println("node count per segment cardinality")
      val perFreqNodeCount = genomeCountPerSegmentWithSegmentLength.groupBy(_._1).mapValues(v => v.size).toList.sortBy(_._1)
      perFreqNodeCount.map { case (x, y) => pw.println(x + "\t" + y) }

      pw.println()
      pw.println("sequence per segment cardinality")
      val perFreqLen = genomeCountPerSegmentWithSegmentLength.groupBy(_._1).mapValues(v => v.map(_._2).sum).toList.sortBy(_._1)
      perFreqLen.map { case (x, y) => pw.println(x + "\t" + y) }

      finish(pw)
      //      pw.close

    }
  }
}