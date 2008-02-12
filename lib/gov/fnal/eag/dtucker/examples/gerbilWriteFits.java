package gov.fnal.eag.dtucker.examples;

import junit.textui.TestRunner;


import nom.tam.image.*;
import nom.tam.util.*;
import nom.tam.fits.*;

/*
 * based on nom.tam.fits.test.ImageTester.java 
 */

public class gerbilWriteFits {
	
    public static void main(String[] args) throws Exception {
		
		Fits f = new Fits();
		
		double[][] dimg = new double[40][40];
		for (int i=0; i<40; i++) {
			for (int j=0; j<40; j++) {
				dimg[i][j] = (double)(i+j);
			}
		}
						
		// Make HDU
		f.addHDU(Fits.makeHDU(dimg));
		
		// Write a FITS file.		
		BufferedFile bf = new BufferedFile("image1.fits", "rw");
		f.write(bf);
		bf.flush();
		bf.close();
		
		System.out.println("finished");

    }
    
}

