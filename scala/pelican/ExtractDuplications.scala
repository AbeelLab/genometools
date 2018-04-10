package pelican

import atk.util.Tool
import java.io.PrintWriter
import java.io.File
import abeel.genometools.Main

object ExtractDuplications extends Main {

  case class Config(val input: File = null, val output: File = null)

  /*
   * /gsap/assembly_repository needs to be mapped to x:
   */
 override def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar pelican.jar extract-duplications") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Input file with ale identifiers and pilon file locations, i.e. ale.id.txt")
      opt[File]('o', "output") action { (x, c) => c.copy(output = x) } text ("Output file. Default input parent directory + duplications.txt")

    } 

    parser.parse(args, Config()).map { config =>
      extract(config);
    }
  }
  def extract(config: Config) {

    val files = tLines(config.input).map(s => s.split("\t")(1).replaceAll("pilon.vcf", "pilon.log"))
    println(files)
    val mapping = files.map(f => {
      println("Processing: " + f)
      f -> tLines(f).filter(_.startsWith("Large collapsed region"))
    })
    //    val oFolder=new File(folder+"/duplications/")
    //    if(!oFolder.exists())
    //      oFolder.mkdirs()
    val pw = new PrintWriter(if(config.output!=null)config.output else new File(config.input.getAbsoluteFile().getParent()+"/duplications.txt"))
    pw.println(generatorInfo)
    pw.println(mapping.map(f => f._1 + "\t" + f._2.mkString(",")).mkString("\n"))
    pw.close
  }

}