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

WARNING: THIS TOOL HAS HARDCODED MAGIC VALUES!!!

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
      val pw= new NixWriter(config.outputFile,config)
      val cc = new Config(config.inputFile, null)

      val gfa = tLines(cc.inputFile)
      //      val pw = new PrintWriter(cc.outputFile)
      val map = gfa.map(line => line.split("\t"))

      
      println("Read in memory")
      val grouped = map.groupBy { x => x(0) }
      val headers = grouped("H")
      pw.println(headers)
      pw.println(headers(1))

      val genomeCount = headers(1)(1).split(";").size
      pw.println("Genome count: " + genomeCount)

//      val h37rvlen = 4411532.0

      val segments = gfa.filter(_(0) == 'S').map(line => {
        val arr = line.split("\t")
        assume(arr(4).startsWith("ORI"))
        //      print(".")
        
        new Segment(arr(1).toInt, arr(2).length(), arr(4).drop(6).replaceAll(".fna", "").replaceAll(".fasta", "").split(";").toList.map(v=>Dictionary.get(v)), MutableList.empty, MutableList.empty)
      })

      pw.println(0)
      pw.println("Mapping identifiers...")
      //map ID to segment
      val nodeMap = segments.map(segment => {
        segment.idx -> segment
      }).toMap

      pw.println("")
      pw.println("Link statistics")
      pw.println("---------------")
      

      /* parse all edges -- has side effect to update segments */
      val links = gfa.filter(_(0) == 'L').map(line => {
        val arr = line.split("\t")
        val from = arr(1).toInt
        val to = arr(3).toInt

        nodeMap(from).outgoing += to
        nodeMap(to).incoming += from

        from.toInt -> to.toInt

      })
      /* Incoming edge distribution */
      
      val inDegree=links.groupBy(_._2).mapValues(_.size)
            
      Histogram.plot(inDegree.toList.map(_._2+0.0),new HistogramConfig(nobin=true,x="In-degree per segment",y="Frequency",outputPrefix=config.outputFile+".indegree",domainStart=0,domainEnd=genomeCount))
      
      
      val histIn=inDegree.map(_._2).groupBy(identity).mapValues(_.size).toList.sortBy(_._1)
      pw.println("In degree distribution: "+histIn.mkString("\n"))
      
      val outDegree=links.groupBy(_._1).mapValues(_.size)
      Histogram.plot(outDegree.toList.map(_._2+0.0),new HistogramConfig(nobin=true,x="Out-degree per segment",y="Frequency",outputPrefix=config.outputFile+".outdegree",domainStart=0,domainEnd=genomeCount))
      val histOut=outDegree.map(_._2).groupBy(identity).mapValues(_.size).toList.sortBy(_._1)
      pw.println("Out degree distribution: "+histOut.mkString("\n"))
      
      pw.println("")
      pw.println("Node statistics")
      pw.println("---------------")
      
      
      
      /* Head nodes */
      val headNodes = segments.filter(seg => seg.incoming.size == 0)
      pw.println("Root nodes: " + headNodes.size)
      pw.println(headNodes.mkString("\n"))

      val tailNodes = segments.filter(seg => seg.outgoing.size == 0)
      pw.println("Tail nodes: " + tailNodes.size)
      pw.println(tailNodes.mkString("\n"))

      nfP.setMaximumFractionDigits(0)

      val ss = grouped("S")
      val ll = grouped("L")
      pw.println("genomes=\t" + genomeCount)
      pw.println("segments=\t" + ss.size)
      pw.println("links=\t" + links.size)
      val sumLen = (ss.map(seg => seg(2).size)).sum
      pw.println("total segment size=\t" + sumLen )
      val genomeCountPerSegmentWithSegmentLength = (ss.map(seg => seg(4).split(";").size -> seg(2).size))
      pw.println("private segment count=\t" + genomeCountPerSegmentWithSegmentLength.filter(_._1 == 1).size)
      pw.println("private segment sum size=\t" + genomeCountPerSegmentWithSegmentLength.filter(_._1 == 1).map(_._2).sum)
      pw.println("private SNP count=\t" + genomeCountPerSegmentWithSegmentLength.filter(_._1 == 1).map(_._2).filter(_ == 1).sum)
      pw.println("shared segment count=\t" + genomeCountPerSegmentWithSegmentLength.filter(f => f._1 > 1 && f._1 < genomeCount).size)
      val fraction = 0.90
      pw.println("shared segment sum size=\t" + genomeCountPerSegmentWithSegmentLength.filter(f => f._1 > 1 && f._1 < (genomeCount * fraction)).map(_._2).sum)

      val hc=new HistogramConfig(nobin=true,x="Genomes per segment",y="Frequency",outputPrefix=config.outputFile+".genomesPerSegment",domainStart=0,domainEnd=genomeCount)
      
      Histogram.plot(genomeCountPerSegmentWithSegmentLength.map(_._1+0.0),hc)

      pw.println("core segment count=\t" + genomeCountPerSegmentWithSegmentLength.filter(f => f._1 == genomeCount).size)
      val coreLen = genomeCountPerSegmentWithSegmentLength.filter(f => f._1 >= (genomeCount * fraction)).map(_._2).sum
      pw.println("core segment sum size=\t" + coreLen )

      finish(pw)
      //      pw.close

    }
  }
}