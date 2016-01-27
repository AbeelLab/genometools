package abeel.genometools

import abeel.genometools.bam2fraglendistr.Bam2FragmentlenDistribution
import abeel.genometools.bam2readnames.Bam2ReadNames
import abeel.genometools.bam2tdf.ConvertBAM2TDF
import abeel.genometools.gbk2gff.GBK2GFF
import abeel.genometools.gff2gtf.GFF2GTF
import abeel.genometools.inject.InjectColumns
import abeel.genometools.reducevcf.ReduceVCF
import abeel.genometools.vcfstats.VCFStatistics
import atk.util.Tool
import abeel.genometools.bamstats.Bamstats
import abeel.genometools.phy2maf.Phy2Maf
import abeel.genometools.sort.MFA
import abeel.genometools.vcf.ConservedRegions
import abeel.genometools.vcf.VCF2gcWindows

trait Main extends Tool {
  def main(args: Array[String]) {}
}


object GenomeToolsConsole extends Tool{

  
  
  def getDeclaredFields(cc: AnyRef) ={
    val m=(Map[String, Any]() /: cc.getClass.getDeclaredFields) { (a, f) =>
      f.setAccessible(true)
      a + (f.getName -> f.get(cc))
    }
    m.toList.sortBy(_._1)
  }
  
  override val version="""
    2015/06/19   Added phy2maf
    2015/07/01   Added vcf2conserved
    2016/01/25   Added vcf2gc
    """
  
  val instructions: Map[String, Main] = Map(

    "bam2fraglendistr" -> Bam2FragmentlenDistribution,

    "bam2readnames" -> Bam2ReadNames,
    "bam2tdf" -> ConvertBAM2TDF,
    "bamstats" ->Bamstats,
    "gbk2gff" -> GBK2GFF,
    "gff2gtf" -> GFF2GTF,
    "inject" -> InjectColumns,
    "reducevcf" -> ReduceVCF,
    "vcfstats" -> VCFStatistics,
    "phy2maf" -> Phy2Maf,
    "sort_mfa" -> MFA,
    "vcf2conserved" -> ConservedRegions,
    "vcf2gc"->VCF2gcWindows
    )
  def main(args: Array[String]): Unit = {

    if (args.length == 0 || !instructions.contains(args(0))) {

      listInstructions
    } else {
      val name = args(0)
      val obj: Main = instructions(name)
      obj.main(args.drop(1))

    }

  }

  def listInstructions() {
    println("Usage:java -jar genometools.jar <instruction> [options...]")
    println("Instructions:")
    println(instructions.toList.sortBy(_._1).map(f => String.format("    %1$-20s",f._1) + f._2.description).mkString("\n"))

  }

}