/*
 * Copyright (c) 2007-2010 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General private License (LGPL), Version 2.1 which
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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package abeel.genometools.bam2tdf;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;


import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

/**
 *   TODO -- normalize option
 *
 -n           Normalize the count by the total number of reads. This option
 multiplies each count by (1,000,000 / total # of reads). It
 is useful when comparing multiple chip-seq experiments when
 the absolute coverage depth is not important.
 */

/**
 * @author jrobinso
 * @author Thomas Abeel
 */
class CoverageCounter {

	private String alignmentFile;
	private Preprocessor consumer;
	private float[] buffer;
	// private int windowSize = 1;
	// TODO -- make mapping qulaity a parameter
	// private int minMappingQuality = 0;
	// FIXME What's this supposed to do?
	// private int strandOption = -1;
	// private int extFactor;
	private int totalCount = 0;

	private SAMSequenceDictionary genome;
	private double mappingQuality=-1;

	CoverageCounter(String alignmentFile, Preprocessor consumer, File wigFile, SAMSequenceDictionary genome2,double mappingQuality) {
		this.alignmentFile = alignmentFile;
		this.consumer = consumer;
		this.genome = genome2;
		buffer = new float[4];
		this.mappingQuality=mappingQuality;
	}

	private boolean passFilter(SAMRecord alignment) {
		return !alignment.getReadUnmappedFlag() && !alignment.getDuplicateReadFlag() && alignment.getMappingQuality()>mappingQuality;
	}

	void parse() throws IOException, URISyntaxException {

	
		SAMFileReader sfr = new SAMFileReader(new File(alignmentFile));

		String lastChr = null;
		ReadCounter counter = null;

//		for (SAMSequenceRecord e : genome.getSequences()) {
			SAMRecordIterator it = sfr.iterator();//sfr.queryOverlapping(e.getSequenceName(), 1, e.getSequenceLength());
			while (it.hasNext()) {
				SAMRecord alignment = it.next();
				if (passFilter(alignment)) {

					totalCount++;

					String alignmentChr = alignment.getReferenceName();//e.getSequenceName();

					if (alignmentChr.equals(lastChr)) {
						if (counter != null) {
							counter.closeBucketsBefore(alignment.getAlignmentStart());
						}
					} else {
						if (counter != null) {
							counter.closeBucketsBefore(Integer.MAX_VALUE);
						}
						counter = new ReadCounter(alignmentChr);
						lastChr = alignmentChr;
					}

					AlignmentBlock[] blocks = alignment.getAlignmentBlocks().toArray(new AlignmentBlock[0]);
					if (blocks != null) {
						for (AlignmentBlock block : blocks) {

							//TODO is this correct?
							int adjustedStart = block.getReferenceStart();// block.getStart();
							int adjustedEnd = block.getReferenceStart() + block.getLength();// block.getEnd();
						
							count(adjustedStart,adjustedEnd,counter,alignment);

							
						}
					} else {
						int adjustedStart = alignment.getAlignmentStart();
						int adjustedEnd = alignment.getAlignmentEnd();
					
						count(adjustedStart,adjustedEnd,counter,alignment);
//						for (int pos = adjustedStart; pos < adjustedEnd; pos++) {
//							if (alignment.getFirstOfPairFlag()) {
//								if (!alignment.getReadNegativeStrandFlag())
//									counter.incrementCount(pos,ReadType.FIRSTREADFORWARDMAP);
//								else
//									counter.incrementCount(pos,ReadType.FIRSTREADREVERSEMAP);
//							}else{
//								if (!alignment.getReadNegativeStrandFlag())
//									counter.incrementCount(pos,ReadType.SECONDREADFORWARDMAP);
//								else
//									counter.incrementCount(pos,ReadType.SECONDREADREVERSEMAP);
//							}
//						}
					}
				}

				

			}
			it.close();
//		}

		if (counter != null) {
			counter.closeBucketsBefore(Integer.MAX_VALUE);
		}

		consumer.setAttribute("totalCount", String.valueOf(totalCount));
		

	}

	

	private void count(int adjustedStart, int adjustedEnd, ReadCounter counter, SAMRecord alignment) {
		
		
		for (int pos = adjustedStart; pos < adjustedEnd; pos++) {
			if (!alignment.getReadPairedFlag()||alignment.getFirstOfPairFlag()) {
				if (!alignment.getReadNegativeStrandFlag())
					counter.incrementCount(pos,ReadType.FIRSTREADFORWARDMAP);
				else
					counter.incrementCount(pos,ReadType.FIRSTREADREVERSEMAP);
			}else{
				if (!alignment.getReadNegativeStrandFlag())
					counter.incrementCount(pos,ReadType.SECONDREADFORWARDMAP);
				else
					counter.incrementCount(pos,ReadType.SECONDREADREVERSEMAP);
			}
		}
		
	}



	class ReadCounter {

		String chr;
		TreeMap<Integer, PositionCounter> counts = new TreeMap<Integer, PositionCounter>();

		ReadCounter(String chr) {
			this.chr = chr;
		}

		void incrementCount(int position, ReadType rt) {
			Integer bucket = position;
			if (!counts.containsKey(bucket)) {
				counts.put(bucket, new PositionCounter());
			}
			counts.get(bucket).increment(rt);
		}

	
		void closeBucketsBefore(int position) {
			List<Integer> bucketsToClose = new ArrayList<Integer>();

			Integer bucket = position;
			for (Map.Entry<Integer, PositionCounter> entry : counts.entrySet()) {
				if (entry.getKey() < bucket) {

					// Divide total count by window size. This is the average
					// count per
					// base over the window, so 30x coverage remains 30x
					// irrespective of window size.
					int bucketStartPosition = entry.getKey();
					int bucketEndPosition = bucketStartPosition + 1;
					if (genome != null) {
						SAMSequenceRecord chromosome = genome.getSequence(chr);
						if (chromosome != null) {
							bucketEndPosition = Math.min(bucketEndPosition, chromosome.getSequenceLength());
						}
					}
					int bucketSize = bucketEndPosition - bucketStartPosition;

					for (int i = 0; i < ReadType.values().length; i++) {
						buffer[i] = ((float) entry.getValue().value(ReadType.values()[i])) / bucketSize;
					}
					

					consumer.addData(chr, bucketStartPosition, bucketEndPosition, buffer, null);

					bucketsToClose.add(entry.getKey());
				}
			}

			for (Integer key : bucketsToClose) {
				counts.remove(key);
			}
			

		}
	}

	private class PositionCounter {

		private int[] counterBuffer = new int[ReadType.values().length];

		void increment(ReadType rt) {
			counterBuffer[rt.ordinal()]++;

		}

		int value(ReadType rt) {
			return counterBuffer[rt.ordinal()];
		}

		
	}

}
