package abeel.genometools.gfa

import abeel.genometools.Main
import java.io.File
import atk.io.NixWriter
import scala.collection.mutable.MutableList
import atk.tools.Histogram.HistogramConfig
import atk.tools.Histogram
import scala.collection.mutable.HashMap
import abeel.genometools.gfa.FastSegment

object FastGFAStatistics extends Main {

  override val description = """Tool to calculate basic statistics from a GFA file.

WARNING: THIS TOOL HAS HARDCODED MAGIC VALUES!!!

This tool is still in development and is not for general use.


"""

  case class Config(val inputFile: File = null, val outputFile: File = null)

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar gfa-fast-statistics") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input GFA formatted file.")
      opt[File]('o', "output") action { (x, c) => c.copy(outputFile = x) } text ("Output file containing statistics.")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)
      val pw = new NixWriter(config.outputFile, config)

      val segments = new HashMap[Int, FastSegment]
      val header = new MutableList[String]
      //      val links = new MutableList[(Int, Int)]
      var linkCount = 0L

      println("reading segments")
      val it = tLinesIterator(config.inputFile)
      while (it.hasNext) {
        progress(1000000)
        val l = it.next()
        //        println(l)

        l(0) match {
          case 'H' => header += l
          case 'S' =>
            val arr = l.split("\t")
            //            assume(arr(4).startsWith("ORI"))
            //      print(".")
            val s = new FastSegment(arr(1).toInt, arr(2).length())
            segments += s.idx -> s

          case 'L' =>
          case 'P' =>
          //            val from = arr(1).toInt
          //            val to = arr(3).toInt
          //
          //            segments(from).outgoing += to
          //            segments(to).incoming += from
          //
          //            links += (from.toInt -> to.toInt)

          case _ => println("Weird line: " + l)
        }
      }
      println("reading links")
      val it2 = tLinesIterator(config.inputFile)
      while (it2.hasNext) {
        progress(1000000)
        val l = it2.next()
        //        println(l)

        l(0) match {
          //          case 'H' => header += l
          //          case 'S' =>
          ////            assume(arr(4).startsWith("ORI"))
          //            //      print(".")
          //            val s = new FastSegment(arr(1).toInt, arr(2).length(), MutableList.empty, MutableList.empty)
          //            segments += s.idx -> s

          case 'L' =>
            val arr = l.split("\t")
            val from = arr(1).toInt
            val to = arr(3).toInt

            segments(from).outgoingCount += 1
            segments(to).incomingCount += 1
            linkCount += 1
          //            links += (from.toInt -> to.toInt)

          case _ =>
        }
      }

      println("Read in memory")
      //      val grouped = map.groupBy { x => x(0) }
      //      val headers = grouped("H")
      //      pw.println(headers)
      //      pw.println(headers(1))

      val genomeCount = 84
      pw.println("")
      pw.println("# Link statistics")
      pw.println("#----------------")

      /* Incoming edge distribution */

      println("in degree")
      val degrees = (segments.map { seg =>( (seg._2.incomingCount + 0.0) ,(seg._2.outgoingCount + 0.0) )}).toList

      //      Histogram.plot(inDegree, new HistogramConfig(nobin = true, x = "In-degree per segment", y = "Frequency", outputPrefix = config.outputFile + ".indegree", domainStart = 0, domainEnd = genomeCount))

      time {

        val histIn = degrees.map(_._1).groupBy(_.toInt).mapValues(_.size).toList.sortBy(_._1)
        pw.println("#In degree distribution: \n" + histIn.map(p => "I\t" + p._1 + "\t" + p._2).mkString("\n"))

      }
      println("out degree")
      time {

        //      Histogram.plot(outDegree, new HistogramConfig(nobin = true, x = "Out-degree per segment", y = "Frequency", outputPrefix = config.outputFile + ".outdegree", domainStart = 0, domainEnd = genomeCount))
        val histOut = degrees.map(_._2).groupBy(_.toInt).mapValues(_.size).toList.sortBy(_._1)
        pw.println("#Out degree distribution: \n" + histOut.map(p => "O\t" + p._1 + "\t" + p._2).mkString("\n"))
      }

      println("Node statistics")
      pw.println("")
      pw.println("# Node statistics")
      pw.println("#----------------")

      /* Head nodes */
      val headNodes = segments.values.filter(seg => seg.incomingCount == 0)
      pw.println("#Root nodes: " + headNodes.size)
      pw.println("#\t" + headNodes.mkString("\n"))

      val tailNodes = segments.values.filter(seg => seg.outgoingCount == 0)
      pw.println("#Tail nodes: " + tailNodes.size)
      pw.println("#\t" + tailNodes.mkString("\n"))

      nfP.setMaximumFractionDigits(0)

      //      val ss = grouped("S")
      //      val ll = grouped("L")
      pw.println("S\tcountGenomes\t" + genomeCount)
      pw.println("S\tcountSegments\t" + segments.size)
      pw.println("S\tcountLinks\t" + linkCount)
      val sumLen = (segments.values.map(seg => seg.sequenceLen)).sum
      pw.println("S\tcumulativeSegmentSize\t" + sumLen)

      //      val genomeCountPerSegmentWithSegmentLength = (segments.values.map(seg => seg.incoming.size -> seg.sequenceLen)).toList
      //      pw.println("private segment count=\t" + genomeCountPerSegmentWithSegmentLength.filter(_._1 == 1).size)
      //      pw.println("private segment sum size=\t" + genomeCountPerSegmentWithSegmentLength.filter(_._1 == 1).map(_._2).sum)
      //      pw.println("private SNP count=\t" + genomeCountPerSegmentWithSegmentLength.filter(_._1 == 1).map(_._2).filter(_ == 1).sum)
      //      pw.println("shared segment count=\t" + genomeCountPerSegmentWithSegmentLength.filter(f => f._1 > 1 && f._1 < genomeCount).size)
      //      val fraction = 0.90
      //      pw.println("shared segment sum size=\t" + genomeCountPerSegmentWithSegmentLength.filter(f => f._1 > 1 && f._1 < (genomeCount * fraction)).map(_._2).sum)

      //      val hc = new HistogramConfig(nobin = true, x = "Genomes per segment", y = "Frequency", outputPrefix = config.outputFile + ".genomesPerSegment", domainStart = 0, domainEnd = genomeCount)
      //
      //      Histogram.plot(genomeCountPerSegmentWithSegmentLength.map(_._1 + 0.0), hc)
      //
      //      pw.println("core segment count=\t" + genomeCountPerSegmentWithSegmentLength.filter(f => f._1 == genomeCount).size)
      //      val coreLen = genomeCountPerSegmentWithSegmentLength.filter(f => f._1 >= (genomeCount * fraction)).map(_._2).sum
      //      pw.println("core segment sum size=\t" + coreLen)

      finish(pw)
      //      pw.close

    }
  }
}