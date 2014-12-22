package abeel.genometools.bamstats

import java.util.Properties
import java.io.File
import atk.util.Tool
import net.sf.samtools.SAMFileReader
import scala.collection.JavaConversions._
import be.abeel.util.FrequencyMap
import be.abeel.util.FrequencyMapUtils
import abeel.genometools.Main
import java.io.PrintWriter

object Bamstats extends Main {

  case class Config(val inputFile: File = null, val outputFile: File = null)

  override val description = "Tool to calculate basic statistics from a BAM file: fragment length distribution, read lenght distribution, read counts."

  override val version = """
    2014/12/22       Initial version
   """

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar bamstats") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input file")
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

    val fm = new FrequencyMap
    val fmReads=new FrequencyMap
    val filtered = sam.iterator().filter(f => f.getFirstOfPairFlag() && !f.getMateUnmappedFlag() && f.getReferenceIndex() == f.getMateReferenceIndex() && f.getMappingQuality() > 0)
    while (filtered.hasNext) {
      val sr = filtered.next
      val s = sr.getAlignmentStart()
      val e = sr.getMateAlignmentStart()
      val diff = math.abs(e - s) + sr.getReadLength()
      fm.count(diff)
      fmReads.count(sr.getReadLength())
    }
    
    val pw=new PrintWriter(config.outputFile+".txt")
    pw.println(generatorInfo)
    pw.println("fraglen.average="+fm.average())
    pw.println("fraglen.median="+fm.median())
    pw.println("fraglen.mode="+fm.mode())
    pw.println("fraglen.count="+fm.totalCount())
    pw.println("fraglen.std="+fm.std())
    pw.println("read.average="+fmReads.average())
    pw.println("read.median="+fmReads.median())
    pw.println("read.mode="+fmReads.mode())
    pw.println("read.count="+fmReads.totalCount())
    pw.println("read.std="+fmReads.std())
    pw.close
    FrequencyMapUtils.plot(fm, config.outputFile.toString()+".png", false, 0, 2000)

    sam.close()

  }

}