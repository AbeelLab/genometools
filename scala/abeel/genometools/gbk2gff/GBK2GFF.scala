package abeel.genometools.gbk2gff;

import java.io.FileOutputStream;
import scala.collection.JavaConversions._
import java.io.IOException;

import net.sf.jannot.Entry;
import net.sf.jannot.EntrySet;
import net.sf.jannot.exception.ReadFailedException;
import net.sf.jannot.parser.FastaParser;
import net.sf.jannot.parser.GFF3Parser;
import net.sf.jannot.parser.Parser;
import net.sf.jannot.source.FileSource;
/**
 * 
 * @author tabeel
 *
 */
object GBK2GFF {

	

	/**
	 * @param args
	 * @throws IOException
	 * @throws ReadFailedException
	 */
	def main(args:Array[String]){
		if (args.length != 2) {
			System.out.println("##----------------------------------------------");
			System.out.println("## GBK2GFF.java");
			System.out.println("## ");
			System.out.println("## Tool to convert a 'gbk' file to a fasta file for");
			System.out.println("## the sequence and a GFF3 file for the annotations.");
			System.out.println("## ");
			System.out.println("## Program will create the output files with the fasta and gff  ");
			System.out.println("## extensions next to the input file.");
			System.out.println("## ");
			System.out.println("## Make sure you have sufficient disk space in");
			System.out.println("## that folder and that you can write there.");
			System.out.println("## ");
			System.out.println("## The program will conclude with a message");
			System.out.println("## that indicates the run was successful.");
			System.out.println("## ");
			System.out.println("## By Thomas Abeel (tabeel@broadinstitute.org)");
			System.out.println("##----------------------------------------------");
			System.out.println("## ");
			System.out.println("## Usage java -jar gbk2gff.jar <input gbk>");
			System.out.println("## ");
		}
		System.out.println("gbk2gff");

		

		val e = new FileSource(args(0)).read();
		val gff = new FileOutputStream(args(0) + ".gff");
		val fasta = new FileOutputStream(args(0) + ".fa");
		for (entry <- e) {
			System.out.println(entry);
			System.out.println(entry.sequence().size());
			Parser.GFF3.write(gff, entry);
			new FastaParser().write(fasta, entry);
		}
		gff.close();
		fasta.close();

	}

}
