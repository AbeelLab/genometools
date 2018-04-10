package pelican

import atk.util.Tool
import be.abeel.util.CountMap
import atk.compbio.vcf.VCFLine
import be.abeel.util.FrequencyMap
import java.io.File
import scala.collection.JavaConversions._
import java.io.PrintWriter
import atk.io.PatternFileFilter
import abeel.genometools.Main

object AggregateAmbiguous extends Main {
  case class Config(val input: File = null, val outputPrefix: String = "ambiguity")
  override def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("java -jar pelican.jar ambiguous") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Folder with VCF files to scan") //, { v: String => config.spacerFile = v })
      opt[String]('o', "output") action { (x, c) => c.copy(outputPrefix = x) } text ("Output prefix")

    }
    parser.parse(args, Config()) map { config =>

      process(config)

    }

  }
  val defaultPattern = ".*.annotated.vcf"
  def process(config: Config) {

    //    val folder = "v:/TB-ARC/KRITH_extended/"
    //    val gs = tLines(folder + "feb18/gnumbers.included.txt")

//    val folder = "v:/TB-ARC/MALI/"
//    val gs = tLines(folder + "gnumbers.included.txt")

    //    val folder = "v:/TB-ARC/ALLAND/"
    //    val gs = tLines(folder + "feb24_30xcontam/gnumbers.included.txt")
    //    
    //    val folder = "v:/TB-ARC/ZHANG/"
    //    val gs = tLines(folder + "gnumbers.included.txt")

    val pw = new PrintWriter(config.outputPrefix + ".summary.txt")
    val pw2 = new PrintWriter(config.outputPrefix + ".perSample.txt")
    val order = List("Amb", "Amb;LowCov", "Del", "Del;Amb", "Del;Amb;LowCov", "Del;LowCov", "LowCov", "PASS")
    pw.println(generatorInfo)
    //    pw.println("# List of files: ")
    //    pw.println("# " + gs.mkString("\n# "))
    val fmArray = Array.ofDim[CountMap[String]](4411709 + 1)
    val actualFm = new CountMap[String]()
    fmArray(0) = new CountMap[String]()
    pw.println("$$\t" + order.mkString("\t") + "\tREFCALL")
    pw2.println("$$\t" + order.mkString("\t"))
    val gs=config.input.listFiles(new PatternFileFilter(defaultPattern))
    for (g <- gs){
     
      val number = g.toString;
       println("Processing: " + g)
      val vcf = config.input+"/" + number + ".annotated.vcf"
      if (new File(vcf).exists) {
        val sizeMap:CountMap[String]=new CountMap[String]
        for(line<-tLinesIterator(vcf).map(new VCFLine(_).filter)){
          sizeMap.count(line)
        }
        
//        val lines = tLines(vcf).toList
//        val sizeMap = lines.map(new VCFLine(_)).groupBy(_.filter).mapValues(_.size)
        
        
        pw2.println(number + "\t" + order.map(sizeMap.getOrElse(_, 0)).mkString("\t"))
        for(line<-tLinesIterator(vcf)){
          val l = new VCFLine(line)
          if (fmArray(l.pos) == null)
            fmArray(l.pos) = new CountMap[String]()
          //          val ft = l.filter.split(";")
          //          ft.map(f => {
          fmArray(0).count(l.filter)
          fmArray(l.pos).count(l.filter)
          //          })
          //          actualFm.count(ft)

        }

        println("# progress: " + fmArray(0))
      } else {
        pw.println("# File missing: " + vcf)
      }

    }
    var count = 0
    for (i <- 1 until fmArray.length) {
      if (fmArray(i) != null) {
        val values = order.map(f => fmArray(i).getOrElse(f, 0).toString.toInt)

        pw.println(i + "\t" + order.map(f => fmArray(i).getOrElse(f, 0)).mkString("\t") + "\t" + (gs.size - values.sum))
        count += 1
      }
    }
    pw.println("# " + count + " positions")
    pw.println("# summary " + fmArray(0))
    val summarizePerType = fmArray(0).keySet().map(key => {
      val arr = fmArray.zipWithIndex.drop(1).filter(p => p._1 != null && p._1.keySet().contains(key))
      val pwx = new PrintWriter(config.outputPrefix + "." + key.replace(';', '_') + ".txt")
      pwx.println(generatorInfo)
      pwx.println("#####")
      pwx.println("## Summary for " + key)
      pwx.println("#####")
      pwx.println(arr.map(f => {
        f._2 + "\t" + f._1.getOrElse(key, -1)
      }).mkString("\n"))
      pwx.close
      key -> arr.size

    })

    pw.println("# positions with particular key")
    pw.println("# " + summarizePerType.map(f => f._1 + "\t" + f._2).mkString("\n# "))
    pw.close
    pw2.close
  }

}