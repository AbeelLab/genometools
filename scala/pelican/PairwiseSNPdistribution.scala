package pelican

import java.io.File
import atk.util.Tool
import be.abeel.util.FrequencyMap
import be.abeel.util.FrequencyMapUtils
import scala.collection.JavaConversions._
import abeel.genometools.Main

/**
 *
 * Program uses input from
 *
 * ~cdesjard/plscripts/fasta_tools/fasta2snpcounts.pl KRITH_extended/feb4.fasta  > KRITH_extended/snpmatrix.txt
 */

object PairwiseSNPdistribution extends Main {
  
  override val version="""
    pre 2015: 	development for TB-ARC
    2015/6/8:	Improved visual output
    
    """
  case class Config(input: File = new File("snp.matrix.txt"), output: File = new File("snp.matrix.distribution.png"))
override  def main(args: Array[String]): Unit = {
 val parser = new scopt.OptionParser[Config]("java -jar pelican.jar snpmatrix-distribution [options]") {
      opt[File]('i', "input")  action { (x, c) => c.copy(input = x) } text ("Input file from SNP matrix program. Default: snp.matrix.txt")
      opt[File]('o', "output") action { (x, c) => c.copy(output = x) } text ("Output file. Default: snp.matrix.distribution.png")
      }

    parser.parse(args, Config()).map { config =>
     run(config)
//      run(new File("v:/TB-ARC/KRITH_extended/snpmatrix.txt"))
   
    }
   

  }

  private def run(config: Config) {

    val map = tMap(tLines(config.input).drop(1))

    val lists = map.values.toList.map(f => (f.split("\t").toList.map(g => g.toInt)).sortBy(identity).drop(1)).flatten

    val binned = lists.map(_.toInt)

    val fm = new FrequencyMap()
    binned.map(fm.count(_))
//    String title, List<FrequencyMap> list, String file,
//			boolean countNormalization, int lower, int upper, String[] labels,
//			String xAxis, String yAxis
    
    FMPlot.plot(null, List(fm), config.output.toString(),false, 0,0,Array("SNPS"),"Pairwise SNP distance","Number of pairs")
    fm.truncate(0, 100)
    FMPlot.plot(null, List(fm), config.output.toString()+"_100",false, 0,0,Array("SNPS"),"Pairwise SNP distance","Number of pairs")
//    FMPlot.plot(fm, config.output.toString()+"_100")
    
  }

}