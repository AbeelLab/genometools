package abeel.genometools

import abeel.genometools.bam2fraglendistr.Bam2FragmentlenDistribution
import abeel.genometools.bam2readnames.Bam2ReadNames
import abeel.genometools.bam2tdf.ConvertBAM2TDF
import abeel.genometools.gbk2gff.GBK2GFF
import abeel.genometools.gff2gtf.GFF2GTF
import abeel.genometools.inject.InjectColumns
import abeel.genometools.reducevcf.ReduceVCF
import abeel.genometools.vcfstats.VCFStatistics
import abeel.genometools.gfa.GFAStatistics

import atk.util.Tool
import abeel.genometools.bamstats.Bamstats
import abeel.genometools.phy2maf.Phy2Maf
import abeel.genometools.sort.MFA
import abeel.genometools.vcf.ConservedRegions
import abeel.genometools.vcf.VCF2gcWindows
import abeel.genometools.vcf.VCF2mutationMatrix
import abeel.genometools.faq.Faq2Kmer
import abeel.genometools.faq.FaqStats
import abeel.genometools.bam.Bam2Kmer
import abeel.genometools.kmer.ReduceKmer
import abeel.genometools.nwk.Nwk2Nodes
import abeel.genometools.nwk.Tree2List
import abeel.genometools.kmer.KmerIntersection
import abeel.genometools.kmer.Kmer2AP
import abeel.genometools.kmer.GroupAP
import abeel.genometools.bam.Bam2GC
import abeel.genometools.faq.Faq2GC
import abeel.genometools.tdf.TDF2GCbias
import abeel.genometools.ena.PrepareENADownload
import net.sf.jannot.mafix.MAFIndex
import abeel.genometools.maf.MAFIndex
import abeel.genometools.vcf.Vcf2MFA


trait Main extends Tool {
  def main(args: Array[String]) {}
  
  
}

object GenomeToolsConsole extends Tool {

  def getDeclaredFields(cc: AnyRef) = {
    val m = (Map[String, Any]() /: cc.getClass.getDeclaredFields) { (a, f) =>
      f.setAccessible(true)
      a + (f.getName -> f.get(cc))
    }
    m.toList.sortBy(_._1)
  }

  override val version = """
    2015/06/19   Added phy2maf
    2015/07/01   Added vcf2conserved
    2016/01/25   Added vcf2gc
    2016/05/23   Added vcf2matrix -- does not work yet
    2016/09/02   Added gfa-statistics
    2016/09/12   Added faq2kmer and faqstats
    2016/09/27   Added bam2kmer, reducekmer, nwk2list, nwk2nodes, intersect-kmer
    2016/09/29   Added bam2gc
    2016/10/05   Added tdf2gcbias
    2016/10/26   Added mafix
    2017/01/13   Added vcf2mfa
    """

  val instructions: Map[String, Main] = Map(

    "bam2fraglendistr" -> Bam2FragmentlenDistribution,

    "bam2readnames" -> Bam2ReadNames,
    "bam2tdf" -> ConvertBAM2TDF,
    "bamstats" -> Bamstats,
    "bam2kmer" -> Bam2Kmer,
    "bam2gc" ->Bam2GC,
    "ena-download"->PrepareENADownload,
    "faq2gc" -> Faq2GC,
    "faq2kmer" -> Faq2Kmer,
    "faqstats" -> FaqStats,
    "gbk2gff" -> GBK2GFF,
    "gff2gtf" -> GFF2GTF,
    "groupprofile"->GroupAP,
    "inject" -> InjectColumns,
    "intersect-kmer" ->KmerIntersection,
    "kmer2matrix"->Kmer2AP,
    "mafix" ->MAFIndex,
    "nwk2list" -> Tree2List,
    "nwk2nodes" -> Nwk2Nodes,
    "reducekmer" -> ReduceKmer,
    "reducevcf" -> ReduceVCF,
    "phy2maf" -> Phy2Maf,
    "sort_mfa" -> MFA,
    "tdf2gcbias" -> TDF2GCbias,
    "vcf2conserved" -> ConservedRegions,
    "vcf2gc" -> VCF2gcWindows,
    "vcf2matrix" -> VCF2mutationMatrix,
    "vcfstats" -> VCFStatistics,
    "vcf2mfa"->Vcf2MFA,
    "gfa-statistics" -> GFAStatistics)
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
    println(instructions.toList.sortBy(_._1).map(f => String.format("    %1$-20s", f._1) + f._2.description).mkString("\n"))

  }

}