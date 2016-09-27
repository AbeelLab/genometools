package abeel.genometools.kmer

import java.io.File
import net.sf.samtools.SAMFileReader
import atk.compbio.DNAString
import scala.collection.JavaConversions._
import java.io.PrintWriter
import abeel.genometools.Main
import net.sf.samtools.SAMFileReader.ValidationStringency
import atk.util.TimeInterval
import atk.compbio.DNAHash

object ReduceKmer extends Main {
  case class Config(val inputFile: File = null, val outputFile: File = null, val count: Int = 5)

  override val description = "Tool to reduce kmer file size by filtering by count."

  override val version = """
    2016/09/27       Initial version included in genometools
   """

  override def main(args: Array[String]) {
 
    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar reducekmer") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input BAM file. ")
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")
      opt[Int]('c', "count") action { (x, c) => c.copy(count= x) } text ("Minimum count to keep, default = " + new Config().count)

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)

      processFile(config)

    }
  }

  private def processFile(config: Config) {

    val pw = new PrintWriter(config.outputFile)
    pw.println(generatorInfo(config))
    
    var counter = 1
    val startTime = System.currentTimeMillis()
    var discard = 0   
    
    for(line<-tLinesIterator(config.inputFile)){
      val arr=line.split("\t")
      val (kmer,count)=arr(0)->arr(1).toInt
      
      
      if (counter % 100000 == 0) {
        val interval = System.currentTimeMillis() - startTime
        println("Processing: " + counter + "\t" + new TimeInterval(interval) + "\t" + nf.format(counter*1000L / (interval + .1)) + " kmers/s")
      }
      counter += 1
 
      
      
      if(count>=config.count){
        pw.println(kmer+"\t"+count)
      }else
        discard += 1
      
    }
    pw.println("# Processed "+counter+" kmers")
    pw.println("# Discarded "+discard+" kmers")
    pw.println("# Done")
    pw.close
    

  }
}