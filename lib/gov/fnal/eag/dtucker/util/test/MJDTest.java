/*
 * Created on Dec 27, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
package gov.fnal.eag.dtucker.util.test;

import gov.fnal.eag.dtucker.util.MJD;
import junit.framework.TestCase;


/**
 * @author dtucker
 *
 */
public class MJDTest extends TestCase {

    public void testGreg2MJD() throws Exception {
        
        double assertToler = 0.0001;
        double[] assertValue1 = new double[10]; 
        double[] assertValue2 = new double[10]; 
        
        // assert values from astrotools tstampToMJD
        assertValue1[0] =    90.641157407;
        assertValue2[0] =   363.641157407;
        assertValue1[1] =  9221.641157407;
        assertValue2[1] =  9495.641157407;
        assertValue1[2] = 18352.641157407;
        assertValue2[2] = 18625.641157407;
        assertValue1[3] = 27483.641157407;
        assertValue2[3] = 27756.641157407;
        assertValue1[4] = 36614.641157407;
        assertValue2[4] = 36887.641157407;
        assertValue1[5] = 45745.641157407;
        assertValue2[5] = 46019.641157407;
        assertValue1[6] = 54877.641157407;
        assertValue2[6] = 55150.641157407;
        assertValue1[7] = 64008.641157407;
        assertValue2[7] = 64281.641157407;
        assertValue1[8] = 73139.641157407;
        assertValue2[8] = 73412.641157407;
        assertValue1[9] = 82270.641157407;
        assertValue2[9] = 82544.641157407;
        
        int YYYY0 = 1859;
        int MM2   = 2;
        int MM11  = 11;
        int DD    = 15;
        int hh    = 15;
        int mm    = 23;
        int ss    = 16;
        
        int YYYY;
        double mjd;
        MJD MJD = new MJD();

        for (int i=0; i<10; i++) {
            YYYY = YYYY0 + 25*i;
            mjd = MJD.greg2MJD( new int[] {YYYY, MM2, DD, hh, mm, ss } );
            System.out.println("MJD for " + YYYY + "-" + MM2 + "-" + DD + 
                    " " + hh + ":" + mm + ":" + ss + " : " + mjd + "  (" + assertValue1[i] + ")");

            assertTrue( (mjd >= assertValue1[i]-assertToler) && 
                        (mjd <= assertValue1[i]+assertToler) ); 
            
            mjd = MJD.greg2MJD( new int[] {YYYY, MM11, DD, hh, mm, ss } );
            System.out.println("MJD for " + YYYY + "-" + MM11 + "-" + DD + 
                    " " + hh + ":" + mm + ":" + ss + " : " + mjd + "  (" + assertValue2[i] + ")");

            assertTrue( (mjd >= assertValue2[i]-assertToler) && 
                        (mjd <= assertValue2[i]+assertToler) ); 

        }
        
    }

    public void testMjd2Greg() throws Exception {

        double[] mjd1 = new double[10];
        double[] mjd2 = new double[10];
        mjd1[0] =    90.641157407;
        mjd2[0] =   363.641157407;
        mjd1[1] =  9221.641157407;
        mjd2[1] =  9495.641157407;
        mjd1[2] = 18352.641157407;
        mjd2[2] = 18625.641157407;
        mjd1[3] = 27483.641157407;
        mjd2[3] = 27756.641157407;
        mjd1[4] = 36614.641157407;
        mjd2[4] = 36887.641157407;
        mjd1[5] = 45745.641157407;
        mjd2[5] = 46019.641157407;
        mjd1[6] = 54877.641157407;
        mjd2[6] = 55150.641157407;
        mjd1[7] = 64008.641157407;
        mjd2[7] = 64281.641157407;
        mjd1[8] = 73139.641157407;
        mjd2[8] = 73412.641157407;
        mjd1[9] = 82270.641157407;
        mjd2[9] = 82544.641157407;

        // assert values from astrotools tstampFromMJD
        int [] assertValueYYYY = new int[] { 1859, 1884, 1909, 1934, 1959, 1984, 2009, 2034, 2059, 2084 }; 
        int assertValueMM1 = 2;  //February
        int assertValueMM2 = 11; //November
        int assertValueDD = 15;  //15:23:16 UT
        int assertValuehh = 15; 
        int assertValuemm = 23; 
        int assertValuess = 16; 

        int assertTolerss = 1;

        
        MJD MJD = new MJD();

        for (int i=0; i<10; i++) {
            
            int[] results1 = MJD.mjd2Greg(mjd1[i]);
            System.out.println("Gregorian Date for MJD " + mjd1[i] + " : " + 
                                results1[0] + "-" + results1[1] + "-" + results1[2] + " " + 
                                results1[3] + ":" + results1[4] + ":" + results1[5] + "  (" + 
                                assertValueYYYY[i] + "-" + assertValueMM1 + "-" + assertValueDD + " " +
                                assertValuehh + ":" + assertValuemm + ":" + assertValuess + ")");
            
            assertEquals(results1[0],assertValueYYYY[i]);
            assertEquals(results1[1],assertValueMM1);
            assertEquals(results1[2],assertValueDD);
            assertEquals(results1[3],assertValuehh);
            assertEquals(results1[4],assertValuemm);
            assertTrue( (results1[5] >= assertValuess-assertTolerss) && 
                        (results1[5] <= assertValuess+assertTolerss) ); 
            
            int[] results2 = MJD.mjd2Greg(mjd2[i]);
            System.out.println("Gregorian Date for MJD " + mjd2[i] + " : " + 
                                results2[0] + "-" + results2[1] + "-" + results2[2] + " " + 
                                results2[3] + ":" + results2[4] + ":" + results2[5] + "  (" + 
                                assertValueYYYY[i] + "-" + assertValueMM2 + "-" + assertValueDD + " " +
                                assertValuehh + ":" + assertValuemm + ":" + assertValuess + ")");
            
            assertEquals(results2[0],assertValueYYYY[i]);
            assertEquals(results2[1],assertValueMM2);
            assertEquals(results2[2],assertValueDD);
            assertEquals(results2[3],assertValuehh);
            assertEquals(results2[4],assertValuemm);
            assertTrue( (results2[5] >= assertValuess-assertTolerss) && 
                        (results2[5] <= assertValuess+assertTolerss) );             

        }

    
    
    }

}
