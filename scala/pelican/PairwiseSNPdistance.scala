package pelican

import net.sf.jannot.source.DataSource
import net.sf.jannot.source.FileSource
import java.io.File
import scala.collection.JavaConversions._
import net.sf.jannot.parser.FastaParser
import java.io.PrintWriter
import atk.util.Tool
import atk.util.TimeInterval
import abeel.genometools.Main

object PairwiseSNPdistance  extends Main{

  case class Config(output: File = new File("snp.matrix.txt"), input: File = null)
 override def main(args: Array[String]): Unit = {
    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar snpmatrix [options]") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Multi-fasta file that's also input for creating phylogenetic tree. ")
      opt[File]('o', "output") action { (x, c) => c.copy(output = x) } text ("Output file. Default: snp.matrix.txt")
    }

    parser.parse(args, Config()).map { config =>
      run(config)
    
    }

  }
  def run(config: Config) {
    FastaParser.forceEntries = true
    val es = new FileSource(config.input).read()
    val pw = new PrintWriter(config.output)
    val pw2 = new PrintWriter(config.output+".comparable.tsv")
    val list = es.iterator().toList.map(e => e.getID() -> e.sequence().stringRepresentation())
    println("Parsed input file")
    val sorted = list.sortBy(_._1)
    pw.print("Taxon\t" + sorted.map(_._1).mkString("\t") + "\n")
    pw2.print("Taxon\t" + sorted.map(_._1).mkString("\t") + "\n")
    var counter=1
    val ts=timestamp
    val start=System.currentTimeMillis()
    for (a <- list) {
      println("Processing: "+a._1 +"\t"+counter+"/"+list.size+"\t"+new TimeInterval(System.currentTimeMillis() - start))
      counter+=1
      pw.print(a._1)
      pw2.print(a._1)
      for (b <- list) {
        
        val paired = a._2.zip(b._2).count(f => f._1 != f._2 && f._1 != 'N' && f._2 != 'N')
        val comparable= a._2.zip(b._2).count(f => f._1 != 'N' && f._2 != 'N')
        pw.print("\t" + paired)
        pw2.print("\t"+comparable)

      }
      pw.print("\n")
       pw2.print("\n")
    }
    pw.close
    pw2.close

  }

}