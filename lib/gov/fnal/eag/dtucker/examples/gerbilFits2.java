package gov.fnal.eag.dtucker.examples;

/*
 * Created on Aug 11, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */

import java.io.*;
import java.util.*;
import nom.tam.util.*;
import nom.tam.fits.*;

/**
 * @author dtucker
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
public class gerbilFits2 {

    public static void main(String[] args) throws Exception {

    		String dirName = "/Users/dtucker/eclipse/desworkspace/gov.fnal.eag.dtucker/gov/fnal/eag/dtucker/examples/data";
        //String fileName = dirName+"/alpha_lyr_stis_002.fits";
        String fileName = dirName+"/p041c_stis_001.fits";
        Fits f = new Fits(fileName);
        
        //BasicHDU[] hdus = f.read();    
        //for (int i=0; i<hdus.length; i += 1) {
        //    Header hdr = hdus[i].getHeader();
        //    System.out.println(hdr.getIntValue("BITPIX"));
        //    hdus[i].info();
        //}

        TableHDU hdu = (TableHDU) f.getHDU(1);
        
        float[] wavelengths = (float[]) hdu.getColumn(0);
        float[] values = (float[]) hdu.getColumn(1);

        System.out.println("file size: " + wavelengths.length + " rows.");
        for (int i=0; i<wavelengths.length; i += 1) {
            System.out.println(wavelengths[i] + "\t" + values[i]);
        }
        
    }
    
}
