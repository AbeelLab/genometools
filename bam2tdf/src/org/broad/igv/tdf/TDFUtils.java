/*
 * Copyright (c) 2007-2010 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tdf;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.URISyntaxException;
import java.util.Map;


import net.sf.samtools.seekablestream.SeekableFileStream;

import org.broad.igv.track.WindowFunction;

/**
 * @author jrobinso
 */
public class TDFUtils {

	public static void dumpRootAttributes(String ibfFile) throws URISyntaxException, FileNotFoundException {
		TDFReader reader = TDFReader.getReader(new SeekableFileStream(new File(ibfFile)));
		System.out.println("Track line = " + reader.getTrackLine());
		TDFGroup group = reader.getGroup("/");
		for (Map.Entry<String, String> entries : group.attributes.entrySet()) {
			System.out.println(entries.getKey() + " = " + entries.getValue());
		}

		System.out.println(reader.getTrackLine());
	}

	public static void dumpDatasets(String ibfFile) throws URISyntaxException, FileNotFoundException {
		TDFReader reader = TDFReader.getReader(new SeekableFileStream(new File(ibfFile)));
		System.out.println("DATASETS");
		for (String dsName : reader.getDatasetNames()) {
			System.out.println(dsName);
			TDFDataset ds = reader.getDataset(dsName);

			System.out.println("Attributes");
			for (Map.Entry<String, String> entry : ds.attributes.entrySet()) {
				System.out.println("\t" + entry.getKey() + " = " + entry.getValue());
			}
			System.out.println();

			System.out.println("Tile Positions");
			for (int i = 0; i < ds.nTiles; i++) {
				System.out.print("\t" + ds.tilePositions[i]);
			}
			System.out.println();

		}
	}

	public static void dumpAllTiles(String ibfFile) throws URISyntaxException, FileNotFoundException {
		TDFReader reader = TDFReader.getReader(new SeekableFileStream(new File(ibfFile)));
		System.out.println("DATASETS");
		for (String dsName : reader.getDatasetNames()) {
			System.out.println(dsName);
			TDFDataset ds = reader.getDataset(dsName);

			for (int i = 0; i < ds.nTiles; i++) {
				TDFTile tile = ds.getTile(i);
				if (tile != null) {
					System.out.println("Tile: " + i);
					dumpTileData(reader, tile);
				}
			}
		}
	}

	public static void dumpTile(String ibfFile, String dsName, int tileNumber) throws URISyntaxException, FileNotFoundException {

		TDFReader reader = TDFReader.getReader(new SeekableFileStream(new File(ibfFile)));
		TDFDataset ds = reader.getDataset(dsName);
		TDFTile tile = reader.readTile(ds, tileNumber);
		if (tile == null) {
			System.out.println("Null tile: " + dsName + " [" + tileNumber + "]");
		} else {
			dumpTileData(reader, tile);
		}
	}

	private static void dumpTileData(TDFReader reader, TDFTile tile) {
		int nTracks = reader.getTrackNames().length;
		int nBins = tile.getSize();
		if (nBins > 0) {
			for (int b = 0; b < nBins; b++) {
				System.out.print(tile.getStartPosition(b));
				System.out.print("\t" + tile.getEndPosition(b));
				System.out.print("\t" + tile.getName(b));
				for (int t = 0; t < nTracks; t++) {
					System.out.print("\t" + tile.getValue(t, b));
				}
				System.out.println();
			}
		}
	}

	public static void dumpRange(String ibfFile, String dsName, int startLocation, int endLocation)
			throws URISyntaxException, FileNotFoundException {
		TDFReader reader = TDFReader.getReader(new SeekableFileStream(new File(ibfFile)));
		TDFDataset ds = reader.getDataset(dsName);

		int tileWidth = ds.tileWidth;
		int startTile = startLocation / tileWidth;
		int endTile = endLocation / tileWidth;

		for (int tileNumber = startTile; tileNumber <= endTile; tileNumber++) {
			TDFTile tile = reader.readTile(ds, tileNumber);
			if (tile == null) {
				System.out.println("Null tile: " + dsName + " [" + tileNumber + "]");
			} else {
				int nTracks = reader.getTrackNames().length;
				int nBins = tile.getSize();
				if (nBins > 0) {
					for (int b = 0; b < nBins; b++) {
						int start = tile.getStartPosition(b);
						int end = tile.getEndPosition(b);
						if (start > endLocation) {
							break;
						}
						if (end >= startLocation) {
							System.out.print(tile.getStartPosition(b));
							for (int t = 0; t < nTracks; t++) {
								System.out.print("\t" + tile.getValue(t, b));
							}
							System.out.println();
						}
					}
				}
			}

		}
	}

	/*
	 * magic number (4 bytes) version index position index size (bytes) # of
	 * window functions [window functions] track type (string) track line
	 * (string) # of tracks [track names]
	 */
	public static void dumpSummary(String ibfFile) throws URISyntaxException, FileNotFoundException {

		TDFReader reader = TDFReader.getReader(new SeekableFileStream(new File(ibfFile)));

		System.out.println("Version: " + reader.getVersion());
		System.out.println("Window Functions");
		for (WindowFunction wf : reader.getWindowFunctions()) {
			System.out.println("\t" + wf.toString());
		}

		System.out.println("Tracks");
		String[] trackNames = reader.getTrackNames();
		for (String trackName : trackNames) {
			System.out.println(trackName);
		}
		System.out.println();

		System.out.println("DATASETS");
		for (String dsName : reader.getDatasetNames()) {
			System.out.println(dsName);
			TDFDataset ds = reader.getDataset(dsName);

			System.out.println("Attributes");
			for (Map.Entry<String, String> entry : ds.attributes.entrySet()) {
				System.out.println("\t" + entry.getKey() + " = " + entry.getValue());
			}
			System.out.println();

			System.out.println("Tiles");

			int nTracks = trackNames.length;
			int tracksToShow = Math.min(4, nTracks);

			for (int i = 0; i < ds.nTiles; i++) {
				TDFTile tile = reader.readTile(ds, i);
				if (tile != null) {
					System.out.print("  " + i);
					/*
					 * int nBins = tile.getSize(); int binsToShow = Math.min(4,
					 * nBins); for (int b = 0; b < binsToShow; b++) {
					 * System.out.print(tile.getStartPosition(b)); for (int t =
					 * 0; t < tracksToShow; t++) { float value =
					 * tile.getValue(0, b); if (!Float.isNaN(value)) {
					 * System.out.print("\t" + tile.getValue(t, b)); } }
					 * System.out.println(); }
					 */
				}
			}
			System.out.println();
			System.out.println();
		}

		System.out.println("GROUPS");
		for (String name : reader.getGroupNames()) {
			System.out.println(name);
			TDFGroup group = reader.getGroup(name);

			System.out.println("Attributes");
			for (Map.Entry<String, String> entry : group.attributes.entrySet()) {
				System.out.println("\t" + entry.getKey() + " = " + entry.getValue());
			}
			System.out.println();

		}
	}

	/*
	 * chr1:241356465-241356657 chr1:241358198-241358223
	 * chr1:241359291-241359329
	 * 
	 * chr4:119691730-119691768 chr4:119692843-119692868
	 * chr4:119694419-119694611
	 */
}
