package abeel.genometools.kmer

import atk.util.Tool
import atk.compbio.DNAHash
import abeel.genometools.Main
import java.io.File
import java.io.PrintWriter
import atk.io.NixWriter
import java.time.LocalDateTime
import scala.collection.mutable.BitSet

class FixedSizeBitSet(nbits: Int) extends BitSet {

  override def toString(): String = {
    val buffer = new StringBuilder(nbits);
    for (i <- 0 until nbits)
      buffer.append(if (this(i)) '1' else '0')

    buffer.toString();
  }
}

object Kmer2AP extends Main {

  case class Config(val outputFile: File = null, val input: File=null)

  override val description = "Generate an A/P matrix for a collection of kmer files"

  override val version = """
    2016/09/27       Initial version included in genometools
   """

  override def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar genometools.jar intersect-kmer") {
      //      opt[File]('i', "input") required () action { (x, c) => c.copy(inputFile = x) } text ("Input BAM file. ")
      opt[File]('o', "output") required () action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("File containing a list of kmer files to be included in the matrix.")

    }
    parser.parse(args, Config()) map { config =>

      assume(config.outputFile != null)
      assume(config.input != null)

      processFile(config)

    }

    //    processFile(new Config(new File("internal_12_0.matrix"), ".bam.reducedkmer", "EscCol_2009EL_2050_ver1_gbk;EscCol_2009EL_2071_ver1_gbk;EscCol_2011C_3493_ver1_gbk"))
  }

  def processFile(config: Config) {

    var hashSize = 0
    val indexedStrains = tLines(config.input).zipWithIndex
    println(indexedStrains)
    indexedStrains.map(f=>assume(new File(f._1).exists,"File does not exist: "+f._1))
    
    val map = scala.collection.mutable.Map[Long, BitSet]() //.withDefaultValue(new FixedSizeBitSet(indexedStrains.size))
    for (strain <- indexedStrains) {
      println("processing: " + strain._1 + "\t" + LocalDateTime.now()+"\t"+map.size)
      for (line <- tLinesIterator(strain._1)) {
        val seq = line.split("\t")(0)
        val ds = DNAHash.hash(seq) //new DNAString(seq)
        if (hashSize == 0)
          hashSize = seq.size
        //        map(ds) += strain._2
        if (!map.contains(ds))
          map += ds -> new FixedSizeBitSet(indexedStrains.size)
        map(ds).add(strain._2)
        //        println(map.size)
      }
    }
    val size = indexedStrains.size
    println("Writing output...")
    val pw = new PrintWriter(config.outputFile)
    pw.println(generatorInfo(config))
    pw.println("files=" + indexedStrains.map(_._1).mkString(";"))

    //    val filtered = map.filter(_._2 == size)

    map.map(f => pw.println(DNAHash.unhash(f._1, hashSize) + "\t" + f._2))
    //    pw.println("# " + filtered.size + " shared kmers")
    pw.close
  }
}


