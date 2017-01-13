package abeel.genometools.bam

import java.io.File
import net.sf.samtools.SAMFileReader
import atk.compbio.DNAString
import scala.collection.JavaConversions._
import java.io.PrintWriter
import abeel.genometools.Main
import net.sf.samtools.SAMFileReader.ValidationStringency
import atk.util.TimeInterval
import atk.compbio.DNAHash
import java.util.Formatter.DateTime
import java.time.LocalDateTime

object Bam2Kmer extends Main{
  case class Config(val inputFile: File = null, val outputFile: File = null, val kmer: Int = 31, val random: Double = 0.01)
  override val description = "Tool to generate subset of kmers and statistics from a BAM file."

  override val version = """
    2016/11/04       Initial version included in genometools
   """
  
  override def main(args: Array[String]) {  
  
    val parser = new scopt.OptionParser[Config]("java -jar mypack.jar thekmers") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input BAM file. ")
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")
      opt[Int]('k', "kmer") action { (x, c) => c.copy(kmer = x) } text ("Kmer length, default = " + new Config().kmer)
      opt[Double]('p', "random") action { (x, c) => c.copy(random = x) } text ("Random value of reads you want processed, default = " + new Config().random)

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)
    
      processFile(config)
      }
    }
  
 

    
    private def processFile(config: Config) {
      val sfr = new SAMFileReader(config.inputFile)
      sfr.setValidationStringency(ValidationStringency.LENIENT)
      val map = scala.collection.mutable.Map[Long, Int]().withDefaultValue(0)
      
      var counter = 0 
      var discard = 0 
  	  val startTime = System.currentTimeMillis() 
  	  
  	  var unmapped = 0
  	  
  	  var skipped = 0 // going to count how many was skipped
  	  val r = scala.util.Random // going to randomly choose which to skip
      val genlength = sfr.getFileHeader().getSequenceDictionary().getReferenceLength()
      var readlengths = 0 // will hold the readlengths to calculate average
      
  	  // iterationg through the reads in the bam file
  	
  	  val fileiter = sfr.iterator() // gets the iterative version of the bam file
 
  	  var howmany = 0 // total reads
    
  	  
  	  while (fileiter.hasNext()) { // will only run till stops
        howmany += 1
     
        if (howmany % 10000 == 0) {
          val interval = System.currentTimeMillis() - startTime
          println("Processing: " + howmany + "\t" + new TimeInterval(interval) + "\t" + nf.format(counter*1000L / (interval + .1)) + " reads/s"+"\t"+LocalDateTime.now())
        }
        if (r.nextFloat > config.random) {
          skipped += 1
          fileiter.next() // skips this
          }
        else {
       	  
          val samRecord = fileiter.next() 
          
          if (samRecord.getReadUnmappedFlag() == true){
            unmapped += 1
          } else {
            counter += 1
            
            val sr = samRecord.getReadString 
            
            readlengths += samRecord.getReadLength()
            
            for (kmerr <- sr.sliding(config.kmer)) { 
              try {
                val ds = DNAHash.hash(kmerr) 
                map(ds) += 1
              }  catch {
                case _: Throwable => discard +=1; 
                }
            }
          }
        }
  	  }
      
      counter = counter - discard
      
      //calculating coverage
      val readlengthavg = readlengths / counter
      val coverage = (counter * readlengthavg) / genlength.toFloat
      
      // printing out the results      
      val pw = new PrintWriter(config.outputFile) // opens new file
      
      pw.println("# Processed reads = "+counter) 
      pw.println("# Discarded reads = "+discard) 
      pw.println("# Skipped reads = "+skipped) 
      pw.println("# Total reads = "+howmany) 
      pw.println("# Coverage = "+coverage) 
      map.map { case (x, y) => pw.println(DNAHash.unhash(x,config.kmer) + "\t" + y)}
      pw.close
      println("Finished")
  }
    
}
