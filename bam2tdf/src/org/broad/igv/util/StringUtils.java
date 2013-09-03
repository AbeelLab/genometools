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

package org.broad.igv.util;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Oct 30, 2009
 * Time: 1:51:36 AM
 * To change this template use File | Settings | File Templates.
 */
public class StringUtils {


    public static List<String> breakQuotedString(String string, char splitToken) {

        ArrayList<String> strings = new ArrayList();
        if (string.length() == 0) {
            return strings;
        }

        char[] characters = string.toCharArray();
        char c = characters[0];

        boolean isQuoted = false;
        StringBuffer buff = new StringBuffer(100);
        for (int i = 0; i < characters.length; i++) {
            c = characters[i];
            if (isQuoted) {
                if (c == '"') {
                    isQuoted = false;
                }
                buff.append(c);
            } else if (c == '"') {
                isQuoted = true;
                buff.append(c);
            } else {
                if (c == splitToken) {
                    strings.add(buff.toString().trim());
                    buff.setLength(0);
                } else {
                    buff.append(c);
                }
            }
        }
        if (buff.length() > 0) {
            strings.add(buff.toString().trim());
        }
        return strings;

    }


    public static short genoToShort(String genotype) {
        byte[] bytes = genotype.getBytes();
        return (short) ((bytes[0] & 0xff) << 8 | (bytes[1] & 0xff));
    }
    
    public static String readString(ByteBuffer byteBuffer) throws IOException {
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        byte b = -1;
        while ((b = byteBuffer.get()) != 0) {
            bytes.write(b);
        }
        return new String(bytes.toByteArray());
    }

    public static void main(String[] args) {

        String genotype = "AC";
        short genoShort = genoToShort(genotype);

        char allel1 = (char) ((genoShort >>> 8) & 0xFF);
        char allel2 = (char) ((genoShort >>> 0) & 0xFF);

        System.out.println("Allele1: " + allel1);
        System.out.println("Allele2: " + allel2);

    }
}
