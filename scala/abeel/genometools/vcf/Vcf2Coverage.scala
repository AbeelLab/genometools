package abeel.genometools.vcf

import scala.io.Source
import java.io.PrintWriter
import java.io.File
import atk.util.Tool
import abeel.genometools.vcf.Mutation._
import abeel.genometools.Main
import atk.compbio.vcf.VCFFile
import scala.collection.JavaConversions._
import be.abeel.util.CountMap
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import atk.io.NixWriter

/**
 * This tool reads VCFs from the path file and outputs a three-column coverage file. Only works if there is a coverage field in the VCF file, which is default for Pilon produced VCF files.
 */
object Vcf2Coverage extends Main {

  case class Config(val vcfFile: File = null, val outputFile: File = null, val window: Int = 1000)

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar vcf2coverage") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(vcfFile = x) } text ("VCF file")
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("Output name.")
      opt[Int]('w', "window") action { (x, c) => c.copy(window = x) } text ("Window size")
    }

    parser.parse(args, Config()) map { config =>
      time {

        val vcf = VCFFile(config.vcfFile)
        val map = scala.collection.mutable.HashMap.empty[String, DescriptiveStatistics]

        for (line <- vcf) {
          progress(10000)
          if (line.pass) {
            val key = line.refGenome + "@@" + (line.pos / config.window)
            if (!map.contains(key))
              map += key -> new DescriptiveStatistics
          if(line.info.contains("DP"))
            map(key).addValue(line.info("DP").toInt)
              
          }
          
        }
        val pw=new NixWriter(config.outputFile.toString,config);
        
        map.toList.sortBy(_._1)(naturalOrdering).map(line=>{
          val id=line._1.split("@@")
          pw.println(id(0)+"\t"+(id(1).toInt*config.window)+"\t"+((id(1).toInt  +1)*config.window-1)+"\t"+line._2.getPercentile(50) +"\t"+line._2.getMean+"\t"+line._2.getValues.size)
        })

        pw.close
        

      }
    }

  }

}