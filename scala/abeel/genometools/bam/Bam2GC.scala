package abeel.genometools.bam

import abeel.genometools.Main
import java.io.File
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileReader.ValidationStringency
import atk.compbio.DNAHash
import scala.collection.JavaConversions._
import atk.io.NixWriter

object Bam2GC extends Main {

  case class Config(val inputFile: File = null, val outputFile: File = null)

  override val description = "Tool to calculate GC bias"

  override val version = """
    2016/09/29       Initial version included in genometools
   """

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar bam2gc") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input BAM file. ")
      opt[File]('o', "output") action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
    

      processFile(config)

    }
  }

  private def gc(c: Byte) = c match {

    case 'a' | 'A' | 't' | 'T' =>
      0;
    case 'c' | 'C' | 'g' | 'G' =>
      1;
    case _ => {
      assert(false, "Cannot encode " + c)
      -1;
    }
  }

  private def processFile(config: Config) {

    val sfr = new SAMFileReader(config.inputFile)
    sfr.setValidationStringency(ValidationStringency.LENIENT)

    val map = scala.collection.mutable.Map[Int, Int]().withDefaultValue(0)
    var discard = 0
    var counter = 0
    for (samRecord <- sfr.iterator()) {
      progress(100000)
      counter += 1

      val sr = samRecord.getReadBases
      //      println(samRecord.getBaseQualities.mkString(""))

      val x = samRecord.getBaseQualities.filter { x => x < 20 }
      //      println("--" + x.mkString(" "))
      if (x.size < 5) {
        val gc = sr.filter(c => c == 'c' | c == 'g' | c == 'G' | c == 'C').size
        map((gc.toDouble / sr.size * 100).toInt) += 1
      } else {
        discard += 1
      }

    }

    val pw = if(config.outputFile!=null) new NixWriter(config.outputFile, config) else new NixWriter(config.inputFile+".gcbias.txt",config)
    pw.println("# Processed reads = " + counter)
    pw.println("# Included reads =" + (counter - discard))
    pw.println("# Discarded reads = " + discard)
    map.toList.sortBy(_._1).map { case (x, y) => pw.println(x + "\t" + y) }

    pw.close

  }
}