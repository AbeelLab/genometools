package abeel.genometools.gff2gtf;

import java.io.FileOutputStream
import scala.collection.JavaConversions._
import java.io.IOException
import net.sf.jannot.Entry
import net.sf.jannot.EntrySet
import net.sf.jannot.exception.ReadFailedException
import net.sf.jannot.parser.GTFParser
import net.sf.jannot.source.FileSource;
import abeel.genometools.Main

/**
 * 
 * @author Thomas Abeel
 * 
 */
object GFF2GTF extends Main{

	/**
	 * @param args
	 * @throws IOException
	 * @throws ReadFailedException
	 */
	override def  main(args:Array[String]){
		if (args.length != 2) {
			System.out.println("##----------------------------------------------");
			System.out.println("## GFF2GTF.java");
			System.out.println("## ");
			System.out.println("## Tool to convert a 'gff' file to a gtf file.");
			System.out.println("## ");
			System.out.println("## Program will create the output file with the gff");
			System.out.println("## extension next to the input file.");
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
			System.out.println("## Usage java -jar gff2gtf.jar <input gbk>");
			System.out.println("## ");
		}
		System.out.println("gff2gtf");

		val e = new FileSource(args(0)).read();
		val gtf = new FileOutputStream(args(0) + ".gtf");
	
		for (entry <- e) {
			System.out.println(entry);
			System.out.println(entry.sequence().size());
			new GTFParser().write(gtf, entry);

		}
		gtf.close();
		

	}

}
