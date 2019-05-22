package pelican

import java.io.File
import atk.util.Tool
import java.io.PrintWriter
import atk.io.PatternFileFilter
import java.nio.file.Files
import java.nio.file.StandardCopyOption
import java.nio.file.CopyOption
import atk.io.IOTools
import abeel.genometools.Main

object MergeProjects extends Main {
override def description = """Merge several projects together."""
  /**
   * Input files contains descriptions of all projects.
   *
   * 1) Top-level directory with outputprefix attached.
   *
   * This program only works on 1 level-deep project nesting!!!
   *
   * It will not deal with merging projects that have overlapping GNumbers!
   *
   * Jobs:
   * Merge following files:
   * - drugs
   * - vcflist
   *
   * Copy files:
   * reduced_vcfs/
   *
   *
   * Create job file:
   * -...
   */

  case class Config(val output: File = null, val input: File = null, val drugSynonyms: File = null, val skipVCF: Boolean = false)

 override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar pelican.jar merge-projects") {
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Output directory prefix (should included dated directory)")
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Input file with description of projects.")
      opt[File]('d', "drug") action { (x, c) => c.copy(drugSynonyms = x) } text ("Input file with synonym descriptions of drugs.")
      opt[Unit]("skip-vcf") action { (x, c) => c.copy(skipVCF = true) } text ("Skip copying VCF files.")

    }

    parser.parse(args, Config()).map { config =>

      merge(config)
    }
  }

  private def merge(config: Config) {
    val lines = tLines(config.input)
    config.output.mkdirs()
    /**
     * Merge drugs
     */
    if (config.drugSynonyms == null) {
      println("WARNING: Drugs will not be merged, requires you to supply a synonym file!")
    }
    if (config.drugSynonyms != null) {
      assume(config.drugSynonyms.exists(), "File does not exist: " + config.drugSynonyms)
      val intermediate = config.output + "/drugs.merge.initial"
      val pwDrug = new PrintWriter(intermediate)
      pwDrug.println(generatorInfo)
      pwDrug.println("# This file is auto-created by merging:")

      val drugsMappingPairs = lines.map(f => {
        val file = f + "/drugs.subset"
        val ls = tLines(file)
        pwDrug.println("#\t" + file)
        val drugs = ls.take(1)(0).split("\t").drop(1)
        /* Map[Gnumber->Map[Drug -> Phenotype]]*/
        val map = tMap(ls.drop(1)).mapValues { f => (drugs.zip(f.split("\t")).toMap) }
        (drugs.toList, map)

      })

      val drugs = drugsMappingPairs.map(_._1).flatten.toSet.toList.sorted(naturalOrdering)
      /* Flattened Map[Gnumber->Map[Drug -> Phenotype]]*/
      val gnumberMap = drugsMappingPairs.map(_._2.toList).flatten.toMap
      val gnumbers = gnumberMap.keySet.toList.sorted(naturalOrdering)
      pwDrug.println("# Total Gnumbers:" + gnumbers.size)
      pwDrug.println("$$\t" + drugs.mkString("\t"))
      for (g <- gnumbers) {
        pwDrug.print(g)
        val drugMap = gnumberMap.getOrElse(g, null)
        assume(drugMap != null)
        for (d <- drugs) {
          pwDrug.print("\t" + drugMap.getOrElse(d, "-"))
        }
        pwDrug.println()

      }
      pwDrug.println("# This file is complete")
      pwDrug.close

      /**
       * Cluster drugs as synonyms
       */
      val pre = tLines(intermediate).map(_.split("\t").toList).transpose.map(f => f(0) -> f).toMap

      val px = new PrintWriter(config.output + "/drugs.subset")
      val clustering = tLines(config.drugSynonyms).map(_.trim.split("\t").toList)

      println(clustering.mkString("\n"))

      val aggregated = (List(pre.getOrElse("$$", null)) ++ clustering.map(drugSynonyms => {
        println(drugSynonyms)
        //rows = drugs, columns are samples , subset of pre
        val drugValues = drugSynonyms.map(drugID => pre.getOrElse(drugID, null))
        println(drugValues)

        //rows= samples, colums =drugs
        val tt = drugValues.filterNot(_ == null).transpose
        val o = tt.map(row => {
          val r = row.filterNot(_.equals("-"))
          if (r.size > 0) {
            println(r)
            if (r(0).size > 1)
              r(0)
            else if (r.contains("1"))
              "1"
            else
              "0"

          } else "-"

        })
        o

      })).transpose

      px.println(aggregated.map(_.mkString("\t")).mkString("\n"))
      px.close

    }
    /**
     * Copy annotated reduced VCF files
     */
    if (!config.skipVCF) {
      val ofx = new File(config.output.getAbsoluteFile().getParentFile() + "/reduced_vcfs/")
      ofx.mkdirs()
      lines.map(f => {
        val pipelineDir = new File(new File(f).getParentFile() + "/reduced_vcfs/")
        /* Default to pipeline directory */
        val fs = if (pipelineDir.exists()) {
          println("Using pipeline directory for " + f)
          new File(new File(f).getParentFile() + "/reduced_vcfs/").listFiles(new PatternFileFilter(".*.annotated.vcf")).par

        } else { /* Else search recursively for *.annotated.vcf files */
          println("Recursively searching for " + f)
          val list = IOTools.recurse(new File(f).getParentFile(), new PatternFileFilter(".*.annotated.vcf"))
          /* Uniqueify */
          list.groupBy(_.getName()).mapValues(_(0)).map(_._2).par

          //        list.par
        }

        println("\tCopying files")
        fs.map(g => {
          val of = new File(ofx + "/" + g.getName())
          if (!of.exists() || of.length != g.length)
            Files.copy(g.toPath(), of.toPath())

        })

      })
    }
    /**
     * Merge files
     */
    List("conversion.txt", "ale.id.txt", "manhattan.id.txt", "bass.id.txt").map { key =>
      val pw = new PrintWriter(config.output.getAbsoluteFile().getParentFile() + "/" + key)
      pw.println(generatorInfo)
      pw.println("# This file is auto-created by merging:")
      val mergedLines = lines.map(f => {
        val file = new File(new File(f).getParentFile() + "/" + key)
        if (!file.exists) {
          pw.println("# Missing file:     " + file)
          List.empty[String]
        } else {
          pw.println("# Included file:    " + file)
          tLines(file)
        }
      })

      pw.println(mergedLines.flatten.mkString("\n"))
      pw.println("# This file is complete")
      pw.close
    }

    List("vcflist.txt.subset", "pgg.txt", "pilon.metrics.txt", "contamination.txt", "ambiguity.perSample.txt", "gnumbers.included.txt").map { key =>

      val pw = new PrintWriter(config.output + "/" + key)
      pw.println(generatorInfo)
      pw.println("# This file is auto-created by merging:")
      val mergedLines = lines.map(f => {
        val file = f + "/" + key
        if (new File(file).exists()) {
          pw.println("# Including file: " + file)
          tLines(file)
        } else {
          pw.println("# Missing file: " + file)
          List.empty[String]
        }
      })

      pw.println(mergedLines.flatten.mkString("\n"))
      pw.println("# This file is complete")
      pw.close
    }
    List("lineages.txt", "spoligotypes.txt").map { key =>
      val pw = new PrintWriter(config.output + "/" + config.output.getName() + "." + key)
      pw.println(generatorInfo)
      pw.println("# This file is auto-created by merging:")
      val mergedLines = lines.map(f => {
        val prefix = new File(f).getName()
        val pipelineFile = new File(f + "/" + prefix + key)
        val file = if (pipelineFile.exists)
          pipelineFile
        else {
          println("Recursive search for .*" + key + " file in " + f)
          val matches = IOTools.recurse(new File(f), new PatternFileFilter(".*" + key))
          assert(matches.size == 1, "Unexpected number of files, please clean up!\n\t" + matches.mkString("\n\t"))
          matches(0)
        }

        pw.println("#\t" + file)
        tLines(file)
      })

      pw.println(mergedLines.flatten.mkString("\n"))
      pw.println("# This file is complete")
      pw.close
    }
  }

}