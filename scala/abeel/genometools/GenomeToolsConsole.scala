package abeel.genometools

import abeel.genometools.bam2fraglendistr.Bam2FragmentlenDistribution
import abeel.genometools.bam2readnames.Bam2ReadNames
import abeel.genometools.bam2tdf.ConvertBAM2TDF
import abeel.genometools.gbk2gff.GBK2GFF
import abeel.genometools.gff2gtf.GFF2GTF
import abeel.genometools.inject.InjectColumns
import abeel.genometools.reducevcf.ReduceVCF
import abeel.genometools.vcfstats.VCFStatistics

trait Main extends Object{
  def main(args: Array[String]){}
}
object GenomeToolsConsole {
  def main(args: Array[String]): Unit = {

    val instructions:Map[String,Main] = Map(

      "bam2fraglendistr" -> Bam2FragmentlenDistribution,

      "bam2readnames" -> Bam2ReadNames,
      "bam2tdf" -> ConvertBAM2TDF,
      "gbk2gff" -> GBK2GFF,
      "gff2gtf" -> GFF2GTF,
      "inject" -> InjectColumns,
      "reducevcf" -> ReduceVCF,
      "vcfstats" -> VCFStatistics)

    if (args.length == 0 || !instructions.contains(args(0))) {

      listInstructions
    } else {
      val name = args(0)
      val obj:Main = instructions(name)
      obj.main(args.drop(1))

    }

  }

  def listInstructions() {
    println("Usage:java -jar tcm.jar [instruction] [instruction options...]")
    println("Instructions:")
    
  }

}