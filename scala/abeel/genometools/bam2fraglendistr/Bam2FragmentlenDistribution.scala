package abeel.genometools.bam2fraglendistr

import java.util.Properties
import java.io.File
import atk.util.Tool
import net.sf.samtools.SAMFileReader
import scala.collection.JavaConversions._
import be.abeel.util.FrequencyMap
import be.abeel.util.FrequencyMapUtils


object Bam2FragmentlenDistribution extends Tool {

  case class Config(val inputFile: File = null, val outputFile: File = null)

  def main(args: Array[String]) {

    println("##----------------------------------------------")
    println("## Bam2FragmentlenDistribution.scala")
    println("## ")
    println("## Tool to calculate a fragment length distribution plot from a BAM file")
    println("## ")
    println("## ")
    println("## By Thomas Abeel (tabeel@broadinstitute.org)")
    println("##----------------------------------------------")

    try {
      val prop = new Properties();
      prop.load(Bam2FragmentlenDistribution.getClass().getResourceAsStream("/tool.properties"));
      println("## Program=" + prop.getProperty("program"));
      println("## Version=" + prop.getProperty("version"));
    } catch {
      case e: Exception =>
        System.err.println("Problem while reading version information");
        e.printStackTrace();
    }

    val parser = new scopt.OptionParser[Config]("java -jar bam2fragmentlenDistribution.jar") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input file in which to inject data") //, { v: String => config.spacerFile = v })
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)

      processFile(config)

    }
  }

  private def processFile(config: Config) {

    

    val sam = new SAMFileReader(config.inputFile)

   
    val fm=new FrequencyMap
    
    val filtered = sam.iterator().filter(f => f.getFirstOfPairFlag() && !f.getMateUnmappedFlag() && f.getReferenceIndex()==f.getMateReferenceIndex() && f.getMappingQuality()>0)
    while (filtered.hasNext){
      val sr=filtered.next
      val s=sr.getAlignmentStart()
      val e=sr.getMateAlignmentStart()
      val diff=math.abs(e-s)+sr.getReadLength()
      fm.count(diff)
    }
    FrequencyMapUtils.plot(fm, config.outputFile.toString(),false,0,2000)


    sam.close()

  }

}