package pelican

import atk.util.Tool
import java.io.File
import java.io.PrintWriter
import abeel.genometools.Main

object PilonMetrics extends Main {
case class Config(input: File = new File("ale.id.txt"), output: File = new File("pilon.metrics.txt"))

 override def main(args: Array[String]): Unit = {
 val parser = new scopt.OptionParser[Config]("java -jar pelican.jar pilon-metrics [options]") {
      opt[File]('i', "input")  action { (x, c) => c.copy(input = x) } text ("Input file from Ale. Default: ale.id.txt")
      opt[File]('o', "output") action { (x, c) => c.copy(output = x) } text ("Output file. Default: pilon.metrics.txt")
      }

    parser.parse(args, Config()).map { config =>
     extract(config.input,config.output)
   
    }
     
  }
  
  private def extract(in:File,out:File){
    val map=tMap(tLines(in))
    val metrics=map.mapValues(_.replaceAll("pilon.vcf", "pilon.metrics"))
    val pw=new PrintWriter(out)
    pw.println(generatorInfo)
    pw.println("# fragCoverage")
    for(line<-metrics){
      val values=tMap(tLines(line._2),limitSplit=false)
      pw.println(line._1+"\t"+values.getOrElse("fragCoverage","-"))
      
    }
    
    pw.close
  }

}