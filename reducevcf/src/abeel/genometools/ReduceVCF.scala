package abeel.genometools

import java.io.File
import java.io.PrintWriter
import java.util.Properties

import scala.io.Source

import atk.util.Tool
import be.abeel.util.CountMap

/**
 * Program to whittle down VCF file to variation sites that pass all filters!
 *
 */
object ReduceVCF extends Tool {

  class Config {
    var inputfile: File = null
    var outputfile: File = null
    var keep: Boolean = false;
  }

  def main(args: Array[String]) {

    println("##----------------------------------------------")
    println("## ReduceVCF.scala")
    println("## ")
    println("## Tool to reduce the size of VCF files by removing")
    println("## all matches.")
    println("## ")
    println("## ")
    println("## The program will conclude with a message")
    println("## that indicates the run was successful.")
    println("## ")
    println("## By Thomas Abeel (tabeel@broadinstitute.org)")
    println("##----------------------------------------------")

    try {
      val prop = new Properties();
      prop.load(ReduceVCF.getClass().getResourceAsStream("/tool.properties"));
      println("## Program=" + prop.getProperty("program"));
      println("## Version=" + prop.getProperty("version"));
    } catch {
      case e: Exception =>
        System.err.println("Problem while reading version information");
        e.printStackTrace();
    }

    val config = new Config();
    val parser = new scopt.OptionParser[Unit]("Reducer") {
      opt[File]('i', "input") required() valueName("<file>") text("Input file") foreach({ v: File => config.inputfile = v })
      opt[File]('o', "output") required() valueName("<file>") text("Output file") foreach({ v: File => config.outputfile = v })
      opt[Unit]('k', "keep") text("Keep all non-reference calls, i.e. calls without the PASS flag.") foreach({ _ => config.keep = true })
      // arglist("<file>...", "arglist allows variable number of arguments",
      //   { v: String => config.files = (v :: config.files).reverse })
    }
    if (parser.parse(args)) {
      assume(config.inputfile!=null)
      assume(config.outputfile!=null)
      processFile(config.inputfile, config.outputfile, config.keep)

    } else {
     
      println("Could not interpret command-line arguments, quitting!")
      System.exit(-1)
    }
  }



  def processFile(file: File, outFile: File, keep: Boolean) = {

   
    if (outFile.exists() && outFile.length() > 0) {
      log("File already exists, aborting...")
    } else {
      val pw = new PrintWriter(outFile)
      val summary = new PrintWriter(outFile.toString() + ".log")
      val filterFM = new CountMap[String]
      val typeFM = new CountMap[String]
      var samePass = 0
      var sameFail = 0
      var diffPass = 0
      var diffFail = 0
      var lineCount = 0
 
      val passCM = new CountMap[String]();
      val failCM = new CountMap[String]();

      for (line <-Source.fromFile(file).getLines) {

        if (line.charAt(0) == '#') {
          pw.println(line)
          //          pwMini.println(line)
        } else {
          lineCount += 1
          val vcfLine = new VCFLine(line)

           filterFM.count(vcfLine.filter)

          typeFM.count(vcfLine.variation.strType)

          if (vcfLine.pass) {
            passCM.count(vcfLine.variation.strType)
          } else {
            failCM.count(vcfLine.variation.strType)
          }

          vcfLine.variation match {
            case p: Match =>
              if (vcfLine.pass)
                samePass += 1
              else {
                if (keep)
                  pw.println(line)
                sameFail += 1
              }
            case _ =>

              if (vcfLine.pass) {
                pw.println(line)
                diffPass += 1
              } else {
                if (keep)
                  pw.println(line)
                diffFail += 1
                //              println("FAIL: " + line)
              }
          }
        
        }

      }

      summary.println(List(file.getName(), samePass, sameFail, diffPass, diffFail, lineCount, passCM, failCM).mkString("\t"))
     
      pw.close()
      summary.close
       finish()
    }

  }

}