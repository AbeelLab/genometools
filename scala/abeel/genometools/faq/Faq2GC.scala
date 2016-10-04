package abeel.genometools.faq

import java.util.Properties
import java.io.File
import atk.util.Tool
import net.sf.samtools.SAMFileReader
import scala.collection.JavaConversions._
import be.abeel.util.FrequencyMap
import be.abeel.util.FrequencyMapUtils
import abeel.genometools.Main
import java.util.HashMap
import atk.compbio.fastq.FastQFile
import atk.compbio.DNAString
import atk.compbio.fastq.FastAFile
import java.io.PrintWriter
import atk.io.NixWriter

object Faq2GC extends Main {

  case class Config(val inputFile: File = null, val outputFile: File = null, val fastq: Boolean = false, val window: Int = 1000)

  override val description = "Tool to draw GC statistics from fasta or fastq file."

  override val version = """
    2016/10/04       Initial version included in genometools
   """

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar faq2gc") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input file. By default FASTA formatted. If you have a FASTQ, use the --fq flag")
      opt[Unit]("fq") action { (x, c) => c.copy(fastq = true) } text ("If you have a FASTQ file, use this flag")
      opt[File]('o', "output") action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")
      opt[Int]('w', "window") action { (x, c) => c.copy(window = x) } text ("Window length, default = " + new Config().window)

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      

      processFile(config)

    }
  }

  private def processFile(config: Config) {

    val gcIt = if (config.fastq)
      FastQFile(config.inputFile).map(fr => fr.seq).map(seq => seq.grouped(config.window))
    else
      FastAFile(config.inputFile).map(fr => fr.seq).map(seq => seq.grouped(config.window))

    val map = scala.collection.mutable.Map[Long, Int]().withDefaultValue(0)
    for (sq <- gcIt.flatten) {
      if (sq.size == config.window) {
        val nt=sq.toUpperCase().groupBy(identity).mapValues { _.size }
        val gc=nt.getOrElse('C', 0)+nt.getOrElse('G', 0)
        val at=nt.getOrElse('A', 0)+nt.getOrElse('T', 0)
        if(at+gc == config.window){
          val fract=math.round((gc*100.0)/(at+gc))
          map(fract) +=1
        }else{
          
        }
        
        
      }

    }

     val pw = if(config.outputFile!=null) new NixWriter(config.outputFile, config) else new NixWriter(config.inputFile+".gc",config)
    
    map.toList.sortBy(_._1).map { case (x, y) => pw.println(x + "\t" + y) }

    pw.close

  }

}