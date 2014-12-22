/*
 * Copyright (c) 2011 by The Broad Institute of MIT and Harvard
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */
package abeel.genometools.bam2tdf

import java.io.File
import java.util.Properties
import org.broad.igv.WindowFunction
import java.util.ArrayList
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileReader.ValidationStringency
import atk.util.Tool
import scala.collection.JavaConversions._
import abeel.genometools.Main

/**
 * Program to create tdf files from bam files.
 *
 * @author Thomas Abeel
 *
 */
object ConvertBAM2TDF extends Tool with Main{
  case class Config(val mappingQuality: Int = -1, files: Seq[File] = Seq())

  override def main(args: Array[String]) {

    println("##----------------------------------------------")
    println("## ConvertBAM2TDF.scala")
    println("## ")
    println("## Tool convert a bam file to a TDF file.")
    println("## ")
    println("## ")
    println("## By Thomas Abeel (tabeel@broadinstitute.org)")
    println("##----------------------------------------------")

    try {
      val prop = new Properties();
      prop.load(ConvertBAM2TDF.getClass().getResourceAsStream("/tool.properties"));
      println("## Program=" + prop.getProperty("program"));
      println("## Version=" + prop.getProperty("version"));
    } catch {
      case e: Exception =>
        System.err.println("Problem while reading version information");
        e.printStackTrace();
    }

    val parser = new scopt.OptionParser[Config]("java -jar bam2tdf.jar") {
      opt[Int]('m', "mappingQuality") action { (x, c) => c.copy(mappingQuality = x) } text ("Minimum mapping quality for reads to be include in the coverage calculation.") //, { v: String => config.spacerFile = v })
      arg[File]("<file>...") unbounded () required () action { (x, c) =>
        c.copy(files = c.files :+ x)
      } text ("files to convert to TDF")

    }
    parser.parse(args, Config()) map { config =>

      for (file <- config.files) {
        processFile(file, config)
      }

    }

  }
  def processFile(file: File, config: Config) = {

    if (!new File(file + ".bai").exists() && !new File(file.toString().replaceAll("\\.bam$", ".bai")).exists()) {

      System.err.println("WARNING: Could not find BAI file for " + file);
      System.err.println("\ttdformat needs a BAI file for each BAM file.");
    } else {
      try {
        createFile(file,config);
      } catch {
        case e:Exception=>
        e.printStackTrace();
        System.err.println("ERROR: Failed to create TDF file for " +file);
      }
    }

  }

  private def printUsage() {
    System.out.println("Usage: java -jar bam2tdf.jar <bam file 1> [<bam file 2> ...]");
    System.out.println("\tbam2tdf needs a BAI file for each BAM file.");

  }

  private def createFile(ifile: File, config: Config) = {
    val wfs = new ArrayList[WindowFunction]();
    for (wf <- WindowFunction.values())
      wfs.add(wf);

    SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
    // TDFTools igvTools = new TDFTools();
    // igvTools.doCount(ifile, ifile + ".tdf", wfs);
    val ofile = ifile + ".tdf";

    System.out.println("Computing coverage.  File = " + ifile);
    val sfr = new SAMFileReader(ifile);
    val dict = sfr.getFileHeader().getSequenceDictionary();
    var max = 0;
    for (ssr <- dict.getSequences()) {
      if (ssr.getSequenceLength() > max)
        max = ssr.getSequenceLength();
    }

    var zoom = 0;
    while (max / 2 > 50000) {
      max /= 2;
      zoom += 1;
    }
    System.out.println("Zoom levels needed: " + zoom);

    val maxZoomValue = zoom;

    val p = new Preprocessor(ifile.getName(), new File(ofile), dict, wfs);
    p.count(ifile.toString(), maxZoomValue, config.mappingQuality);
    p.finish();

    System.out.flush();

  }

}
