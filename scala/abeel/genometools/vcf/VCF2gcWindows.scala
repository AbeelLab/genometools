package abeel.genometools.vcf

import atk.util.Tool
import atk.compbio.vcf.VCFFile
import atk.util.TimeInterval
import java.io.File
import java.io.PrintWriter
import abeel.genometools.Main

object VCF2gcWindows extends Main {

  case class Config(val input: File = null, val output: File = null, val window:Int=1000, val incomplete:Boolean=false)

  override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar vcf2gc") {
      opt[File]('i', "input") action { (x, c) => c.copy(input = x) } text ("Input file")
      opt[File]('o', "output") action { (x, c) => c.copy(output = x) } text ("Output file")
      opt[Int]("window") action { (x, c) => c.copy(window= x) } text ("Window size")
      opt[Unit]("incomplete")action { (x, c) => c.copy(incomplete= true) } text ("Include incomplete windows")

    }
    parser.parse(args, Config()) map { config =>

   
      var time = System.currentTimeMillis()
      val file = VCFFile(config.input)
      val linesIn = file.toList
      println("read: " + new TimeInterval(System.currentTimeMillis() - time))
      time = System.currentTimeMillis()
     

      //      pwLog.println("filtering: " + new TimeInterval(System.currentTimeMillis() - time))
      //      pwLog.flush()
      time = System.currentTimeMillis()
      val linesX = linesIn.par.filter(_.pass).filter(l => l.altLength == 1 && l.refLength == 1)
      println("Nfilter: " + new TimeInterval(System.currentTimeMillis() - time))
    
      val lines = linesX.filterNot(l => l.ref.contains("N")).toList
      println(lines.take(10).mkString("\n"))
  
      println("group per chr")
     
      val pw = new PrintWriter(config.output)
      pw.println(generatorInfo)
      pw.println("##gcRate\tavgDepths\tat\tgc\tvalues\tdepth values")
      var complete = 0
      var incomplete = 0
      var verify = 0
      val chrs = lines.groupBy(_.refGenome)
      for ((chr, cLines) <- chrs) {
        println("processing " + chr + "\t" + cLines.size)
     
        val xx = cLines.groupBy(_.zeroPos / config.window)

        for ((id, window) <- xx) {

          if (window.size != config.window) {
            incomplete += 1
            println("WARNING: INCOMPLETE WINDOW!")
          } else {
            complete += 1
            
            
          }

          val at = window.filter(l => l.ref.equals("A") || l.ref.equals("T"))
          val gc = window.filter(l => l.ref.equals("G") || l.ref.equals("C"))
          val atCount = at.toList.size
          val gcCount = gc.toList.size
          if (atCount + gcCount == config.window)
            verify += 1
          val gcRate = gcCount / (atCount.toDouble + gcCount.toDouble)
          val dp = (at ++ gc).map(l => l.arr(7).split(";")(0).split("=")(1).toInt)

          //          assume(atCount + gcCount == xx.size)

          if(atCount + gcCount == config.window || config.incomplete)
        	  pw.println(gcRate + "\t" + (dp.sum / dp.size.toDouble) + "\t" + at.size + "\t" + gc.size + "\t" + (at.size + gc.size) + "\t" + dp.size)
        }
      }
      
      pw.println("# incomplete windows=" + incomplete)
      pw.println("# complete windows=" + complete)
      pw.println("# verify=" + verify)
      pw.close()

    }
  }

}