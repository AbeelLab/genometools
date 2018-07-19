package pelican

import atk.util.Tool
import java.io.PrintWriter
import java.io.File
import abeel.genometools.Main

/**
 * Tool to flag suspect samples based on:
 *
 *  - Pilon fragment coverage
 *  - Ambiguous variant calls
 *  - Contamination
 *  - CPT rate between LongInsertions and SingleSubstitions
 *
 *
 *   Future:
 *
 *   -  SNP distance compared to neighbors
 *   	* Calculate baseline distribution by taking all subsequent pairs in tree and calculating the distance
 *    	between their SNP-distance profile. Normal distribution with mean and variance.
 *      * Calculate location of sample within this distribution,if outside N std -> report
 */

object FlagSuspectSamples extends Main {
  case class Config(val manhattan: String = "1", val threshold: Double = 1, val outputPrefix: String = "", val ambiguous: File = null, val pilonMetrics: File = null, val contamination: File = null, val cptTable: File = null)
 
  override def description = """Identify samples that should be excluded based on their CPT, "+
        "ambiguity rate and large insertion rate."""
  
  override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar pelican.jar drugs") {
      opt[File]("ambiguity") required () action { (x, c) => c.copy(ambiguous = x) } text ("Input file with ambiguity information per sample.")
      opt[File]("pilon-metrics") required () action { (x, c) => c.copy(pilonMetrics = x) } text ("Input file with pilon metrics.")
      opt[File]("contamination") required () action { (x, c) => c.copy(contamination = x) } text ("Input file with contamination information.")
      opt[File]("cpt-table") required () action { (x, c) => c.copy(cptTable = x) } text ("Input file with CPT table.")
      opt[Double]("threshold") action { (x, c) => c.copy(threshold = x) } text ("Score threshold to exclude samples. Default = 1 ")
      opt[String]('o', "output") action { (x, c) => c.copy(outputPrefix = x) } text ("Output prefix ")
      opt[String]('m', " manhattan") required () action { (x, c) => c.copy(manhattan = x) } text ("Manhattan identifiers as a comma separated list")

    }

    parser.parse(args, Config()).map { config =>
      run(config);
    }

  }
  def run(config: Config) {

    val amb = tMap(tLines(config.ambiguous))
    val cov = tMap(tLines(config.pilonMetrics))
    val contam = tMap(tLines(config.contamination))

    val cpt = tMap(tLines(config.cptTable))

    val keys = amb.keys ++ cov.keys ++ contam.keys
    assume(amb.getOrElse("$$", null).split("\t")(0).equals("Amb"))
    assume(amb.getOrElse("$$", null).split("\t")(7).equals("PASS"))
    nf.setMaximumFractionDigits(3)
    val pw = new PrintWriter(config.outputPrefix + "suspicious.all.txt")
    val pw2 = new PrintWriter(config.outputPrefix + "suspicious.removed.txt")
    val pw3 = new PrintWriter(config.outputPrefix + "suspicious.retained.txt")
    val tuples = for (key <- keys.filterNot(_.equals("$$"))) yield {
      println(key)
      if (!amb.contains(key)) {
        pw.println("# WARNING: missing ambiguity rate for key " + key)
        pw2.println("# WARNING: missing ambiguity rate for key " + key)
        pw3.println("# WARNING: missing ambiguity rate for key " + key)
      }
      val ambArr = if (amb.contains(key)) amb.getOrElse(key, null).split("\t").map(_.toDouble) else null
      val ambRate = if (ambArr == null || ambArr(7) < 25) 0 else ambArr(0) / ambArr(7)

      val lisArr = cpt.getOrElse(key, null).split("\t").map(_.toDouble)
      val liRate = if (lisArr(0) < 25) 0 else lisArr(5) / lisArr(0)

      val lsRate = if (lisArr(0) < 25) 0 else lisArr(3) / lisArr(0)

      val covV = cov.getOrElse(key, "0").toDouble.toInt

      val contaminant = contam.getOrElse(key, "-").trim
      var suspicionScore = 0.0
      //      if (ambRate > 0.5)
      suspicionScore += ambRate
      //      if (liRate > 0.5)
      suspicionScore += liRate
      suspicionScore += lsRate
      if (covV < 30)
        suspicionScore += 5
      //      if (contaminant.size > 1)
      //        suspicionScore += 1

      //      pw.println(key + "\t" + nf.format(ambRate) + "\t" + nf.format(liRate) + "\t" + covV + "\t" + contaminant + "\t" + nf.format(suspicionScore))
      (key, ambRate, liRate, lsRate, covV, contaminant, suspicionScore)
    }

    pw.println(generatorInfo)
    pw2.println(generatorInfo)
    pw3.println(generatorInfo)
    pw.println("#")
    pw.println("## Configuration: ")
    pw.println("# ambiguous = " + config.ambiguous)
    pw.println("# contamination = " + config.contamination)
    pw.println("# pilon-metrics = " + config.pilonMetrics)
    pw.println("# cpt-table = " + config.cptTable)
    pw.println("# threshold = " + config.threshold)
    pw.println("################################")
    pw.println("# Samples marked for exclusion")
    pw.println("################################")
    pw.println("## List available in: " + config.outputPrefix + "suspicious.removed.txt")
    pw.println("##G\tAmbiguity rate\tLargeInsertion rate\tLargeSubstitution rate\tfragCoverage\tcontaminati\tscore")
    tuples.filter(_._7 >= config.threshold).map(f => printTuple(pw, f))
    tuples.filter(_._7 >= config.threshold).map(f => printTuple(pw2, f))
    pw.println("################################")
    pw.println("# Samples included")
    pw.println("################################")
    pw.println("##G\tAmbiguity rate\tLargeInsertion rate\tLargeSubstitution rate\tfragCoverage\tcontamination\tscore")
    tuples.filter(_._7 < config.threshold).map(f => printTuple(pw, f))
    tuples.filter(_._7 < config.threshold).map(f => printTuple(pw3, f))

    def printTuple(pwx: PrintWriter, tuple: (String, Double, Double, Double, Int, String, Double)) {
      val (key, ambRate, lsRate, liRate, covV, contaminant, suspicionScore) = tuple
      pwx.println(key + "\t" + key + "\t" + nf.format(ambRate) + "\t" + nf.format(liRate) + "\t" + nf.format(lsRate) + "\t" + covV + "\t" + contaminant + "\t" + nf.format(suspicionScore))
    }

    pw.close
    pw2.close
    pw3.close
  }

}
