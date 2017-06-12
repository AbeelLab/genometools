package abeel.genometools.tdf

import net.sf.samtools.SAMFileReader
import java.io.PrintWriter
import java.io.File
import abeel.genometools.Main
import java.time.LocalDateTime
import atk.compbio.DNAHash
import net.sf.samtools.SAMFileReader.ValidationStringency
import atk.util.TimeInterval
import scala.collection.JavaConversions._
import net.sf.jannot.tdf.TDFDataSource
import net.sf.jannot.source.Locator
import net.sf.jannot.source.DataSourceFactory
import net.sf.jannot.tdf.TDFData
import atk.io.NixWriter
import java.util.logging.Logger
import java.util.logging.Level
import atk.util.LoggingTrait

object TDF2GCbias extends Main {

 
  
  case class Config(val input: File = null, val reference: File = null, val window: Int = 1000, val output: File = null)

  override val description = "Tool to calculate GC bias from coverage file (TDF) and a reference (fasta)."

  override val version = """
    2016/10/05       Initial version included in genometools
   """

  override def main(args: Array[String]) {
 setDebugLevel(Level.WARNING)
    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar tdf2gcbias") {
      val default = new Config()
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Input TDF file. ")
      opt[File]('r', "reference") required () action { (x, c) => c.copy(reference = x) } text ("Input fasta file. ")
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("File where you want the output to be written")
      opt[Int]('w', "window") action { (x, c) => c.copy(window = x) } text ("Window length, default = " + default.window)

    }
    parser.parse(args, Config()) map { config =>

      //      assume(config.inputFile != null)
      //      assume(config.outputFile != null)

      processFile(config)

    }
  }

  private def processFile(config: Config) {

    val ref = DataSourceFactory.create(new Locator(config.reference)).read();
    val tdf = DataSourceFactory.create(new Locator(config.input)).read();

    val pw = new NixWriter(config.output, config)

    pw.println("# GC-content per window of "+config.window+" nt")
    pw.println("#reference\tcoverage")
    for (entry <- ref) {
      pw.println("#ENTRY = "+entry.getID)
      println(entry.getID + "\t" + entry.iterator().toList)
      val r = entry
      val tEntry = tdf.getEntry(entry.getID)

      val t = (tEntry.get(tEntry.iterator().toList(0))).asInstanceOf[TDFData]
      for (i <- 0 until r.getMaximumLength / config.window) {
        val seq = r.sequence().subsequence(i * config.window, (i + 1) * config.window).stringRepresentation()
        val nt = seq.toUpperCase().groupBy(identity).mapValues { _.size }
        val gc = nt.getOrElse('C', 0) + nt.getOrElse('G', 0)
        val at = nt.getOrElse('A', 0) + nt.getOrElse('T', 0)

        val cov = t.get(i * config.window, (i + 1) * config.window).toList
        val covMap = cov.map { pile => (pile.start() -> pile.end()) -> pile.getTotal }
        val singles=covMap.filter(p=>p._1._1+1==p._1._2).filter(p=>p._1._1>=i * config.window && p._1._1<(i + 1) * config.window)
        val singleSum=singles.map(_._2).sum
//        println(covMap.size+"\t"+singles.size)
//        covMap.filter(p)
        if (at + gc == config.window) {
          val fract = math.round((gc * 100.0) / (at + gc))
          val cc = singleSum / config.window
          pw.println(fract + "\t" + cc)
//        } else {
//          pw.println(fract + "\t" + cc)
//          println("Incomplete window: " + (at + gc) + ", " + singles.size)
        }

      }

    }
    pw.close

  }
}