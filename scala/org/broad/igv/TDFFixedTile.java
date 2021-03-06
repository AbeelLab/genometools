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
package org.broad.igv;

import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * @author jrobinso
 */
public class TDFFixedTile implements TDFTile {
	/* (non-Javadoc)
	 * @see org.broad.igv.tdf.TDFTile#noValues()
	 */
	@Override
	public int noValues() {
		return data.length;
	}
    int tileStart;
    float span;
    int start;
    float[][] data;

    public TDFFixedTile(ByteBuffer byteBuffer, int nSamples) throws IOException {
        this.fill(byteBuffer, nSamples);
    }

    public TDFFixedTile(int tileStart, int start, float span, float[][] data) {
        this.tileStart = tileStart;
        this.span = span;
        this.data = data;
        this.start = start;
    }

    public int getTileStart() {
        return start;
    }

    public int getTileEnd() {
        return getSize() == 0 ? 0 : getEndPosition(getSize() - 1);
    }

    public int getStartPosition(int idx) {
        return start + (int) (idx * span);
    }

    public int getEndPosition(int idx) {
        return start + (int) ((idx + 1) * span);
    }

    public String getName(int idx) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public float getValue(int row, int idx) {
        return data[row][idx];
    }

    public int getSize() {
        return (data == null ? 0 : data[0].length);
    }


    // TODO -- record "type",  extent (longest feature), other stuff

    public void writeTo(BufferedByteWriter fos) throws IOException {

        fos.putNullTerminatedString(TDFTile.Type.fixedStep.toString());
        fos.putInt(getSize());
        fos.putInt(start);
        fos.putFloat(span);
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                fos.putFloat(data[i][j]);
            }
        }

    }

    public void fill(ByteBuffer byteBuffer, int nSamples) throws IOException {

        int nPositions = byteBuffer.getInt();
        start = byteBuffer.getInt();
        span = byteBuffer.getFloat();

        data = new float[nSamples][nPositions];
        for (int sample = 0; sample < nSamples; sample++) {
            data[sample] = new float[nPositions];
            for (int i = 0; i < nPositions; i++) {
                data[sample][i] = byteBuffer.getFloat();
            }
        }

    }
}
