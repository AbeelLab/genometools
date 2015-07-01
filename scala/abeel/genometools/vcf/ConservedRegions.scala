package abeel.genometools.vcf

import java.io.File
import java.io.PrintWriter
import atk.compbio.vcf.VCFLine
import be.abeel.util.FrequencyMap
import be.abeel.util.FrequencyMapUtils
import atk.util.Tool
import abeel.genometools.GenomeToolsConsole
import abeel.genometools.Main

object ConservedRegions extends Main {

  override val description = """
    Tool to identify perfectly conserved regions from a set of VCF files.
    
    Currently only works for single chromosome genomes.
    
    """

  case class Config(val includeSNV: Boolean = false, val reportLen: Int = 100, val genomeLen: Int = 0, val input: File = null, val output: File = new File("conserved_regions.txt"), val flankLSV: Int = 100, val flankSNV: Int = 100)

  override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar genometools vcf2conserved") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("File that contains the list of VCF files to scan")
      opt[File]('o', "output") action { (x, c) => c.copy(output = x) } text ("File where you want the output to be written")
      opt[Int]("lsv-flank") action { (x, c) => c.copy(flankLSV = x) } text ("How many flanking bases are also removed from conserved area for each LSV. (default=100)")
      opt[Int]("snv-flank") action { (x, c) => c.copy(flankSNV = x) } text ("How many flanking bases are also removed from conserved area for each SNV. (default=100). Note that by default SNV's are not taking into account, see 'include-SNV' option.")
      opt[Unit]("include-SNV") action { (x, c) => c.copy(includeSNV = true) } text ("Include SNVs to determine conserved regions.")
      opt[Int]('l', "len") required () action { (x, c) => c.copy(genomeLen = x) } text ("Length of the genome.")
      opt[Int]("min-len") action { (x, c) => c.copy(reportLen = x) } text ("Minimum length of reported regions. (default=100)")

    }
    parser.parse(args, Config()) map { config =>

      group(config)

    }
  }

  private def group(config: Config) {

    val list = tLines(config.input).map(new File(_))

    val pw = new PrintWriter(config.output + ".intermediate")
    pw.println(generatorInfo)
    pw.println("# Tool configuration: ")
    pw.println("# " + GenomeToolsConsole.getDeclaredFields(config).mkString("\n# "))

    val lsvArr = Array.ofDim[Int](config.genomeLen)
    val snpArr = Array.ofDim[Int](config.genomeLen)

    pw.println("# Files included: " + list.size)

    for (file <- list) {
      println("Processing: " + (list.indexOf(file) + 1) + "/" + list.size + " " + file)
      val lines = tLinesIterator(file)

      while (lines.hasNext) {
        val line = lines.next
        val vcf = new VCFLine(line)
        val start = vcf.pos
        val len = vcf.refLength
        if (vcf.passOrUnfiltered) {
          if (vcf.altLength > 1 || vcf.refLength > 1) {
            for (i <- start - config.flankLSV to start - config.flankLSV + len) {
              if (i >= 0 && i < lsvArr.size) {
                lsvArr(i) += 1
              }
            }
          } else {
            for (i <- start - config.flankSNV to start - config.flankSNV + len) {
              if (i >= 0 && i < lsvArr.size) {
                snpArr(i) += 1
              }
            }
          }
        }
      }

    }

    //    pw.println("# SNP bases: " + snpArr.filter(_ > 0).size)
    pw.println("# LSV bases: " + lsvArr.filter(_ > 0).size)
    pw.println("# ")
    pw.println("## SNP-density\tLSV-density")

    val zip = snpArr.zip(lsvArr)
    pw.println(zip.map(f => f._1 + "\t" + f._2).mkString("\n"))

    pw.close

    /**
     * Convert information to actual regions
     */
    val values = tLines(config.output + ".intermediate").map(_.split("\t").map(_.toInt))
    val lsv = values.map(f => f(1) + (if (config.includeSNV) f(0) else 0))

    def encode(s: List[Int]) = {
      val coded = s.foldLeft((0, s(0), List.empty[(Int, Int)]))((t, char) =>
        t match {
          case (i, p, output) => if (p == char) (i + 1, p, output) else (1, char, (i -> p) :: output)

        }) match { case (i, p, output) => (i -> p) :: output }
      coded.reverse
    }

    val rle = encode(lsv)
    val startPositions = rle.foldLeft(List(0))((last, current) =>
      last.head + current._1 :: last).reverse

    val merged = startPositions.zip(rle)
    val filtered = merged.filter(p => p._2._2 == 0 && p._2._1 >= config.reportLen)
    val pwx = new PrintWriter(config.output)
    pwx.println(generatorInfo)
    pwx.println("# Number of regions: " + filtered.size)
    pwx.println("## position\tlength")
    pwx.println(filtered.map(f => f._1 + "\t" + f._2._1).mkString("\n"))
    pwx.close
    println("done")

  }

}