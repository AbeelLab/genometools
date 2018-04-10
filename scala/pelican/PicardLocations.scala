package pelican

import java.io.File
import atk.util.Tool
import java.io.PrintWriter
import atk.io.DirectoryFilter
import abeel.genometools.Main

/**
 * Make sure \\iodine-cifs\seq_picard_aggregation is mapped to K
 *
 * K itself will not be readable, but each G subdirectory will be
 */

object PicardLocations extends Main {
  case class Config(val input: File = new File("manhattan.id.txt"), val output: File = new File("bass.id.txt"), val root: String = "/seq/picard_aggregation/")
  override def main(args: Array[String]): Unit = {

    //    findBASS("k:/", new File("v:/TB-ARC/development/manhattan.id.txt"), new File("v:/TB-ARC/development/bass.test.txt"))

    val parser = new scopt.OptionParser[Config]("java -jar pelican.jar picard") {
      opt[File]('i', "input") action { (x, c) => c.copy(input = x) } text ("Input file from Manhattan identifier extraction")
      opt[File]('o', "output") action { (x, c) => c.copy(output = x) } text ("Output file.")
      opt[String]('r', "root") action { (x, c) => c.copy(root = x) } text ("Root directory to search from, by default should be /seq/picard_aggregation/")

    }

    parser.parse(args, Config()).map { config =>
      findBASS(config.root, config.input, config.output);
    }

  }

  def findBASS(root: String, input: File, output: File) = {
    val gs = tLines(input).map(_.split("\t")).map(_(2).split(";").toList).flatten.filterNot(_.equals("-"))
    println(gs)
    val pw = new PrintWriter(output)
    pw.println(generatorInfo)
    pw.println("# root=" + root)
    pw.println("# input=" + input)
    pw.println("# g\tlibrary_type\tlocation\taggregateLocation\tINITIATIVE	PROJECT	SAMPLE_ALIAS	PAIRED_RUN	FLOWCELL_BARCODE	LANE	LIBRARY_NAME	ALIGNED_BAM	REFERENCE_SEQUENCE	TARGET_INTERVALS	BAIT_INTERVALS	LIBRARY_TYPE	MOLECULAR_BARCODE_SEQUENCE	SAMPLE_LSID	INDIVIDUAL_ALIAS	MOLECULAR_BARCODE_NAME")
    gs.map(g => {
      println("finding: " + g)
      val dir = new File(root + "/" + g)
      if (dir.exists) {
        println("found " + dir)
        println(dir.listFiles().toList)
        val subdir = dir.listFiles(new DirectoryFilter).filter(f => f.listFiles().length > 0 && !f.getName().startsWith("gssr"))
        assume(subdir.size == 1)
        val allversions = subdir(0).listFiles()
        val versions = if (allversions.exists(_.getName().equals("current"))) {
          println("Using current")
          allversions.filter(_.getName().equals("current"))
        } else {
          println("Current not found")
          allversions.sortBy(f => f.getName().substring(1).toInt)
        }

        println("Versions: " + versions.toList)
        val lines = tLines(versions.last + "/analysis_files.txt")
        val agg = new File(versions.last + "/" + versions.last.getParentFile().getName() + ".bam")
        assume(agg.exists(), "AGG does not exist: " + agg)

        val libs = lines.drop(1).map(_.split("\t"))
        libs.map(arr =>
          pw.println(g + "\t" + arr(11) + "\t" + arr(7) + "\t" + agg + "\t" + arr.mkString("\t")))
      } else {
    	  pw.println("# "+g+"\tNo files found in Picard")
      }
    })

    pw.close
  }

}