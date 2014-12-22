package abeel.genometools.inject

import java.util.Properties
import java.io.File
import atk.util.Tool
import java.io.PrintWriter

object InjectColumns extends Tool{

  case class Config(val inputFile: File = null, val outputFile: File = null, files: List[File] = List())

  def main(args: Array[String]) {

    println("##----------------------------------------------")
    println("## InjectColumns.scala")
    println("## ")
    println("## Tool to inject extra columns in a tab-delimited file based on key-value")
    println("## ")
    println("## ")
    println("## The program will conclude with a message")
    println("## that indicates the run was successful.")
    println("## ")
    println("## By Thomas Abeel (tabeel@broadinstitute.org)")
    println("##----------------------------------------------")

    try {
      val prop = new Properties();
      prop.load(InjectColumns.getClass().getResourceAsStream("/tool.properties"));
      println("## Program=" + prop.getProperty("program"));
      println("## Version=" + prop.getProperty("version"));
    } catch {
      case e: Exception =>
        System.err.println("Problem while reading version information");
        e.printStackTrace();
    }

    val parser = new scopt.OptionParser[Config]("java -jar injectcolumns.jar") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input file in which to inject data") //, { v: String => config.spacerFile = v })
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")

      //      opt("o", "output", "<file>", "File where you want the output to be written", { v: String => config.outputFile = v })
      //
      arg[File]("<file>...") unbounded () required () action { (x, c) => c.copy(files = c.files :+ x) } text ("List of files that contains data you want injected. They need to be key\tvalue pairs. These files are applied in the order specified. So later files may inject based on information injected in earlier files.")

      //      arg("<file>", "<file> is the BAM file you like to spoligotype", { v: String => config.inputFile = v })
      // arglist("<file>...", "arglist allows variable number of arguments",
      //   { v: String => config.files = (v :: config.files).reverse })
    }
    parser.parse(args, Config()) map { config =>

      assume(config.inputFile != null)
      assume(config.outputFile != null)
      
      processFile(config)

    }
  }
  
  private def processFile(config:Config){
    
    
    var content=tLines(config.inputFile,false,false).mkString("\n")
    
    config.files.map(f=>{
    	val kv=tMap(tLines(f))
    	kv.map(p=>{
    	  content=content.replaceAll(p._1, p._1+"\t"+p._2)
    	})
    	
      
    })
    
    val pw=new PrintWriter(config.outputFile)
    pw.println(content)
    pw.close
    
    
    
  } 
  
}