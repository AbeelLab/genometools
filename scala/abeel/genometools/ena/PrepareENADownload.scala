package abeel.genometools.ena

import atk.io.URLCache
import java.io.PrintWriter
import atk.util.Tool
import java.io.File
import abeel.genometools.Main

object PrepareENADownload extends Main {

  override val version = """
    Tool to prepare downloads from ENA from a list of identifiers.
    
    2014/11/19: Initial version
    2014/12/02: Added option to organize files for BI specific projects
                Include collaborator identifiers when available for BI projects
    2014/12/09: Fixed sample identifier in the first conversion file column to 
                reflect difference between Broad and other projects 
    2015/10/15: If ENA files are missing try to dump from SRA
    2016/09/12: Initial version in genometools
   """

  case class Config(val input: String = null, val debug: Boolean = false, val output: File = null, val broad: Boolean = false)

  override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar ena-download") {
      opt[String]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Single ENA identifier or filename with list of identifiers.") //, { v: String => config.spacerFile = v })
      opt[File]('o', "output") action { (x, c) => c.copy(output = x) } text ("Output file, by default output is written to console.")
      opt[Unit]("debug") action { (x, c) => c.copy(debug = true) } text ("Show debug output")
      opt[Unit]("broad") action { (x, c) => c.copy(broad = true) } text ("Use Broad Institute organization of projects, with an individual project per sample.")

    }
    parser.parse(args, Config()) map { config =>
      val f = new File(config.input)
      if (f.exists()) {
        prepFromList(tColumn(0, tLines(f)), config.output, config)
      } else
        prepFromList(List(config.input), config.output, config)
    }

  }

  def prepFromList(list: List[String], of: File, config: Config) = {
    val f = new File(config.input)
    val pw = if (of != null) {
      of.getAbsoluteFile().getParentFile().mkdirs()
      new PrintWriter(of+".download.sh")
    } else {
      new PrintWriter(System.out)
    }

    val idMap = if (f.exists) {
      if (tLines(f)(0).split("\t").size > 1) {
        tMap(tLines(f), keyColumn=0,valueColumn=1, limitSplit = false)
      } else {
        tLines(f).zip(tLines(f)).toMap
      }

    } else Map.empty[String, String]

    pw.println(generatorInfo(config))
    pw.println("## Data download instructions for " + list)

    val pwConversion =

      if (of != null) {

        new PrintWriter(of + ".conversion.txt")
      } else {
        new PrintWriter(System.out)
      }

    pwConversion.println(generatorInfo(config))

    val pwx = if (of != null) {

      new PrintWriter(of + ".sourcequery.txt")
    } else {
      new PrintWriter(System.out)
    }
    pwx.println(generatorInfo(config))

    for (l <- list) {

      val url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + l + "&result=read_run&fields=run_accession,library_layout,fastq_ftp,submitted_ftp,experiment_alias,scientific_name,sample_alias&download=text"
      println("URL: " + url)
      val data = URLCache.query(url)

      pwx.println("## URL=" + url)
      pwx.println(data.mkString("\n"))

      //    val samples=tLines(new URL(url).openStream())
      val samples = data.drop(1)
      //      for (sample <- samples) {

      pw.println("# #samples: " + samples.size)
      //    pw.println(samples.mkString("\n"))
      if(samples.size==0){
        
      }
      
      samples.map(s => {
        val arr = s.split("\t")
        if (config.debug)
          println("SAMPLE: "+s)

        val fdr = if (config.broad) l else arr(0)
        pw.print("mkdir " + fdr + "\n")
        pw.print("cd " + fdr + "\n")
        val sampleNameExpand = if (arr(3).trim.size > 0)
          arr(3).split("/").last.replaceAll(".fastq.bz2", "")
        else
          arr(4).split("_").dropRight(2).drop(1).mkString("_")

        val sampleName = if (arr(3).trim.size > 0 && arr(3).split(";").size > 1)
          sampleNameExpand.dropRight(2)
        else sampleNameExpand
        if (config.debug){
          println("L: "+l+"\tIDmap: "+idMap)
        }
        pwConversion.println(fdr + "\t" + arr(0) + "\t" + sampleName.replaceAll("_sequence", "-") + "\t" + idMap.getOrElse(l, "-")+"\t"+arr(4)+";"+arr(5)+";"+arr(6))
        val column = if (arr(2).trim.size == 0) arr(3) else arr(2)
        val files = column.split(";")
        
        files.filterNot(_.size==0).map(f => pw.print("wget -q " + f + "\n"))
        pw.print("cd ..\n")
        
        if (files(0).size == 0) {
          //assume(files(0).size > 0, "There are no FASTQ files available for " + arr(0))
          pw.print("## There are no FASTQ files available for " + arr(0) + ", let's try direct SRA dump\n")
          
          pw.print("../bin/fastq-dump --split-files -O "+fdr+"/"+" "+arr(0)+"\n")
          print("## There are no FASTQ files available for " + arr(0) + "\n")
        }
        
       

      })
      pw.flush()
      pwx.flush()
      pwConversion.flush()
    }
    //      }
    pw.close
    pwx.close
    pwConversion.close

  }

}