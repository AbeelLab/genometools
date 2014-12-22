package abeel.genometools.gff2gtf;

import java.io.FileOutputStream;
import java.io.IOException;

import net.sf.jannot.Entry;
import net.sf.jannot.EntrySet;
import net.sf.jannot.exception.ReadFailedException;
import net.sf.jannot.parser.GTFParser;
import net.sf.jannot.source.FileSource;

/**
 * 
 * @author Thomas Abeel
 * 
 */
public class GFF2GTF {

	/**
	 * @param args
	 * @throws IOException
	 * @throws ReadFailedException
	 */
	public static void main(String[] args) throws ReadFailedException, IOException {
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

		EntrySet e = new FileSource(args[0]).read();
		FileOutputStream gtf = new FileOutputStream(args[0] + ".gtf");
	
		for (Entry entry : e) {
			System.out.println(entry);
			System.out.println(entry.sequence().size());
			new GTFParser().write(gtf, entry);

		}
		gtf.close();
		

	}

}
