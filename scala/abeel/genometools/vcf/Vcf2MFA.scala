package abeel.genometools.vcf

import scala.io.Source
import java.io.PrintWriter
import java.io.File
import atk.util.Tool
import abeel.genometools.vcf.Mutation._
import abeel.genometools.Main

/**
 * This tool reads VCFs from the path file, and returns a SNP phy-file.
 * It does not remove samples with duplicate names.
 */
object Vcf2MFA extends Main {

  case class Config(val vcfPathFile: File = null, val output: File = null, val mq:Int=0)

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar vcf2mfa") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(vcfPathFile = x) } text ("VCF path file")
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Output name.")
      opt[Int]("mq") action {(x,c)=>c.copy(mq=x)} text("Filter on mapping quality < <value> (default = 0, no filtering)")
    }

    

    /* Get all SNP positions */
    def getPositions(fList: List[File]): Map[String, String] = {
      def getPos(file: File): List[(String, String)] = {
        //        val sample = file.getAbsoluteFile.getParentFile
        println("Reading " + file.getParentFile().getName+"/"+file.getName + "...")
        val snpIterator = Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filterNot(_.isConfirmedReference).filter(_.isSNP)
        val snpList=snpIterator.map(_ match {
          case SNP(r, c, a, chr) => (chr + "_" + c, r)
        }).toList
        snpList
      }
      fList.flatMap(getPos(_)).toMap
    }

    //    /* Return new name of samples with 9 characters and a space. */
    //    def truncateName(s: String): String = {
    //      if (s.length > 10) s.substring(s.length - 9, s.length) + " " //Cut s to length 10
    //      else { val res = "          ".substring(0, 10 - s.length); s + res } //Add spaces until length 10
    //    }

    /**
     *  List all VCFs from the path file, get all SNP positions,
     *  and print SNP sequence of each sample.
     */

    parser.parse(args, Config()) map { config =>
      time {
        SNP.mqFilter=config.mq
        val vcfList = tLines(config.vcfPathFile).map(new File(_)).toList
        val refMap = getPositions(vcfList) // Map with ref. positions and bases
        println("Writing mfa-file...")
        val p = new PrintWriter(config.output)
        val pw = new PrintWriter(config.output + ".log")
        pw.println(generatorInfo())
        pw.println("Number of VCF files: " + vcfList.size)
        pw.println("Positions in file: " + refMap.size)
        pw.println("Positions in aligments: ")
        pw.println(refMap.keySet.toList.sorted.mkString("\n"))

        //        p.println(vcfList.size + 1 + " " + refMap.size) //Print total number of sequences (VCF's) + reference (1st sequence) and total number of SNP positions.
        /* Output constructed reference */
        p.println(">reference")
        refMap.keysIterator.toList.sorted.foreach(pos => p.print(refMap(pos)))
        p.println
        
        
        vcfList.foreach { file => // for each file print sequence
          val name = file.getName
          val snpMap = tLines(file).filterNot(_.isConfirmedReference).filter(_.isSNP).map(_ match {
            case SNP(r, c, a, chr) => (chr + "_" + c, a)
          }).toMap
          
          val nonSnpSet = tLines(file).filterNot(_.isConfirmedReference).filter { line =>
            val arr = line.split("\t")
            val id = arr(0) + "_" + arr(1)
            refMap.contains(id)
          }.map {
            line =>
              // line => line.split("\t")(1).toInt
              val arr = line.split("\t")
              val id = arr(0) + "_" + arr(1)
              id
          }.toSet

          val snpSeq = refMap.keysIterator.toList.sorted.map(pos =>
            if (snpMap.contains(pos)) snpMap(pos)
            else if (nonSnpSet.contains(pos)) "N"
            else refMap(pos)).mkString
          p.println(">" + name)
          p.println(snpSeq)
          println(name + ":\t" + snpMap.size + "\tSNPs")
        }
        p.close
        pw.println("Total of " + vcfList.size + " VCFs read.")
        pw.println("Length of SNP sequences: " + refMap.size)
        pw.println("Output: " + args(1))
        pw.close
      }
    }

  }

}