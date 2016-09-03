package abeel.genometools.gfa
import abeel.genometools.Main
import java.io.File
import java.io.PrintWriter
import scala.collection.mutable.MutableList

/**
 * Program to calculate statistics on GFA files
 */

object GFAStatistics extends Main {

  override val description = """Tool to calculate basic statistics from a GFA file.

WARNING: THIS TOOL HAS HARDCODED MAGIC VALUES!!!

This tool is still in development and is not for general use.


"""

  case class Config(val inputFile: File = null, val outputFile: File = null)

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar gfa-statistics") {
            opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input GFA formatted file.")
//            opt[File]('o', "output") action { (x, c) => c.copy(outputFile = x) } text ("Output file containing statistics.")

    }
    parser.parse(args, Config()) map { config =>

            assume(config.inputFile != null)
      //      assume(config.outputFile != null)

      val cc = new Config(config.inputFile, null)

      val gfa = tLines(cc.inputFile)
//      val pw = new PrintWriter(cc.outputFile)
      val map = gfa.map(line => line.split("\t"))

      val grouped = map.groupBy { x => x(0) }
      val headers = grouped("H")
      println(headers)
      println(headers(1))

      val genomeCount = headers(1)(1).split(";").size
      println("Genome count: " + genomeCount)

      val h37rvlen = 4411532.0

      val segments = gfa.filter(_(0) == 'S').map(line => {
        val arr = line.split("\t")
        assume(arr(4).startsWith("ORI"))
        //      print(".")
        new Segment(arr(1).toInt, arr(2), arr(4).drop(6).replaceAll(".fna", "").replaceAll(".fasta", "").split(";").toList, MutableList.empty, MutableList.empty)
      })

      println(0)
      println("Mapping identifiers...")
      //map ID to segment
      val nodeMap = segments.map(segment => {
        segment.idx -> segment
      }).toMap

      println("Parsing links...")

      /* parse all edges -- has side effect to update segments */
      val links = gfa.filter(_(0) == 'L').map(line => {
        val arr = line.split("\t")
        val from = arr(1).toInt
        val to = arr(3).toInt

        nodeMap(from).outgoing += to
        nodeMap(to).incoming += from

        from.toInt -> to.toInt

      })

      val headNodes = segments.filter(seg => seg.incoming.size == 0)
      println("Root nodes: " + headNodes.size)
      println(headNodes.mkString("\n"))

      val tailNodes = segments.filter(seg => seg.outgoing.size == 0)
      println("Tail nodes: " + tailNodes.size)
      println(tailNodes.mkString("\n"))
      
      nfP.setMaximumFractionDigits(0)
      
      val ss = grouped("S")
      val ll = grouped("L")
      println("genomes=\t"+genomeCount)
      println("segments=\t" + segments.size)
      println("links=\t" + links.size)
      val sumLen = (ss.map(seg => seg(2).size)).sum
      println("total segment size=\t" + sumLen + " (" + nfP.format(sumLen / h37rvlen)+")")
      val genomeCountPerSegmentWithSegmentLength = (ss.map(seg => seg(4).split(";").size -> seg(2).size))
      println("private segment count=\t" + genomeCountPerSegmentWithSegmentLength.filter(_._1 == 1).size)
      println("private segment sum size=\t" + genomeCountPerSegmentWithSegmentLength.filter(_._1 == 1).map(_._2).sum)
      println("private SNP count=\t" + genomeCountPerSegmentWithSegmentLength.filter(_._1 == 1).map(_._2).filter(_ == 1).sum)
      println("shared segment count=\t" + genomeCountPerSegmentWithSegmentLength.filter(f => f._1 > 1 && f._1 < genomeCount).size)
      val fraction = 1.0
      println("shared segment sum size=\t" + genomeCountPerSegmentWithSegmentLength.filter(f => f._1 > 1 && f._1 < (genomeCount * fraction)).map(_._2).sum)
      
      
 
      println("core segment count=\t" + genomeCountPerSegmentWithSegmentLength.filter(f => f._1 == genomeCount).size)
      val coreLen = genomeCountPerSegmentWithSegmentLength.filter(f => f._1 >= (genomeCount * fraction)).map(_._2).sum
      println("core segment sum size=\t" + coreLen + " (" + nfP.format(coreLen / h37rvlen)+")")


      
      
//      pw.close

    }
  }
}