package pelican

import atk.util.Tool
import atk.io.PatternFileFilter
import java.io.File
import atk.compbio.vcf.VCFFile
import atk.compbio.vcf.VCFLine
import atk.util.TimeInterval

import atk.util.NaturalOrderComparator
import java.io.PrintWriter
import abeel.genometools.Main

object LSP2genes extends Main {

  class GFFLine(val line: String) {
    lazy val arr = line.split("\t")
    //    assume(arr.length)
    lazy val `type` = arr(2)
    lazy val start = arr(3).toInt
    lazy val end = arr(4).toInt
    lazy val len = end - start + 1
    lazy val notes = arr(8).split(";").toList.map { f => val pair = f.split("="); pair(0).toLowerCase() -> pair(1) }.toMap
    def field(key: String) = notes.getOrElse(key.toLowerCase, null)
    def id() = field("ID")
  }

 override def main(args: Array[String]): Unit = {
    val parser = new scopt.OptionParser[Config]("java -jar pelican.jar variant2gene") {
      opt[String]('o', "output") required() action { (x, c) => c.copy(output = x) } text ("Output prefix")
      opt[File]("vcf") required() action { (x, c) => c.copy(vcfDirectory= x) } text ("Directory with all VCF files. Need to end with '.annotated.vcf'")
      opt[File]("gff") required() action { (x, c) => c.copy(gff= x) } text ("File with annotations in GFF format")
      opt[String]("type") action{(x,c)=>c.copy(typ=x)}text("Type (third column in gff) of feature to use, default=gene")
      opt[String]("id") action{(x,c)=>c.copy(id=x)}text("Identifier to use for display (last GFF column), default=ID")
      
      opt[File]('i',"input")required() action { (x, c) => c.copy(includList= x) } text ("File with list of G-numbers to include")
      
    }

    parser.parse(args, Config()).map { config =>
      lsp2genes(config);
    }
  }

  case class Config(val id:String="ID",gff: File = null, includList: File = null, vcfDirectory: File = null, output: String = null,val typ:String="gene")
  def lsp2genes(config: Config) {
    val genes = (tLines(config.gff).map(new GFFLine(_)).filter(_.`type`.equalsIgnoreCase(config.typ))).par
    println("Genes included in analysis: "+genes.size)
    val includeList = tLines(config.includList)//.map(GNumbers.singleG(_))
    val list = config.vcfDirectory.listFiles(new PatternFileFilter(".*.annotated.vcf")).filter(file => includeList.contains(file.getName)).toList.par
    //    val list=new File("v:/TB-ARC/KRITH_extended/reduced_vcfs/").listFiles(new PatternFileFilter(".*.annotated.vcf")).toList.par
    val t = System.currentTimeMillis()
    val data = list.map { file =>
      val g = file.getName
      println("Processing: " + file)
      val hitsPerLine = VCFFile(file).filter(f => f.pass && !f.alt.equals(".")).map(vcf => {
        val hits = genes.filter(overlap(vcf, _))
        vcf -> hits
      })
//           println("\t"+hitsPerLine.length)
      val triplets = hitsPerLine.toList.map(pair => pair._2.map(f => (f.field(config.id), pair._1, g)).toList)
      triplets
    }
    println("Read all data")
    val flatdata = data.flatten.flatten

    println("Flattened")
    val perRV = flatdata.groupBy(_._1)

    val perRV_type = perRV.mapValues(f => f.groupBy(_._2.variation.toString))

    println("Grouped: ")
    val sortedType = flatdata.groupBy(_._2.variation.toString).keys.toList.sorted(naturalOrdering)
    println("Type list: "+sortedType)
    val sortedRv = perRV_type.keys.toList.sorted(naturalOrdering)
    val pw = new PrintWriter(config.output + "variant2gene.txt")
    pw.println(generatorInfo)
    pw.println("$$\t" + sortedType.mkString("\t") + "\tSingleEvent\tLongEvent\tAnyEvent\t" + sortedType.mkString("\t") + "\tSingleEvent\tLongEvent\tAnyEvent")
    for (rv <- sortedRv) {
      pw.print(rv)

      val sets = sortedType.map { t =>
        val query = perRV_type.getOrElse(rv, null).getOrElse(t, null)
        val set = if (query == null) Set.empty[String] else query.map(_._3).toSet
        pw.print("\t" + (set.size))
        t->set
      }.toSet

//      pw.print("\t" + sets.flatten.size)
      pw.print("\t"+
          sets.filter(x=>x._1.startsWith("Single")).map(_._2).flatten.size+"\t"+  
          sets.filter(x=>x._1.startsWith("Long")).map(_._2).flatten.size+"\t" + 
          sets.map(_._2).flatten.size)

      val setsExpanded = sortedType.map { t =>
        val query = perRV_type.getOrElse(rv, null).getOrElse(t, null)
        val set = if (query == null) Set.empty[String] else query.map(_._3).toSet
        //    	if(set.size>200){
        //    	  println("Oops: "+set.mkString("\n"))
        //    	 System.exit(0) 
        //    	}
        pw.print("\t" + (if (set == null) "-" else set.mkString(",")))
        t->set
      }
      pw.print("\t"+
          setsExpanded.filter(x=>x._1.startsWith("Single")).map(_._2).flatten.mkString(",")+"\t"+  
          setsExpanded.filter(x=>x._1.startsWith("Long")).map(_._2).flatten.mkString(",")+"\t" +
          setsExpanded.map(_._2).flatten.mkString(","))
      pw.println()
    }

    pw.close

    println("Time: " + new TimeInterval(System.currentTimeMillis() - t))

  }

  def overlap(vcf: VCFLine, gff: GFFLine): Boolean = {
    (vcf.pos <= gff.end) && (vcf.pos + vcf.refLength >= gff.start)
  }

}