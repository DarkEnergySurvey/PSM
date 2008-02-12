package gov.fnal.eag.dtucker.util;

/*
 * The class MJD contains menthods for converting between 
 * Gregorian Dates and MJD for dates after 1858 Nov 17 00:00:00 UT.
 * 
 * Ref:  Chapter 7 of Jean Meeus's Astronomical Algorithms, 
 *       2nd Edition
 * 
 * Created on Dec 22, 2005
 *
 */

public class MJD {
    
    /**
     * @author dtucker
     * @param YMDhms
     * Convert from Gregorian calendar date (YYYY MM DD hh mm ss) to MJD, 
     * for dates after 1858 Nov 17 00:00:00 UT.
     * (Conversion accurate to about 1 sec in time.)
     * 
     */

    public double greg2MJD(int[] YMDhms) throws Exception {

        int YYYY = YMDhms[0];
        int MM   = YMDhms[1];
        int DD   = YMDhms[2];
        int hh   = YMDhms[3];
        int mm   = YMDhms[4];
        int ss   = YMDhms[5];
        
        if (MM == 1 || MM == 2) {
            YYYY = YYYY - 1;
            MM = MM + 12;
        }
        
        double dDD = DD + (hh + (mm + ss/60.)/60.) / 24.;
        
        int a = (int) (YYYY / 100.);
        int b = 2 - a + (int) (a / 4.);
        double jd = (int) (365.25*(YYYY+4716)) + (int)(30.6001*(MM+1)) + dDD + b -1524.5;
     
        double mjd = jd-2400000.5;

        if (mjd < 0) {
            throw new Exception();
        }
        
        return mjd;
        
    }
    
    /**
     * @author dtucker
     * @param MJD
     * Convert from MJD to Gregorian calendar date (YYYY MM DD hh mm ss), 
     * for dates after 1858 Nov 17 00:00:00 UT.
     * (Conversion accurate to about 1 sec in time.)
     * 
     */
    public int[] mjd2Greg(double mjd) throws Exception {
        
        if (mjd < 0) {
            throw new Exception();
        }
        
        double jd = mjd+2400000.5;       
        
        double jdp = jd + 0.5;
        int      z = (int) jdp;
        double   f = jdp - z;
        
        int a;
        if (z < 2299161) {
            a = z;
        } else {
            int alpha = (int) ( (z - 1867216.25) / 36524.25 );
            a = z + 1 + alpha - (int) (alpha/4.);
        }
        
        int b = a + 1524;
        int c = (int) ( (b - 122.1) / 365.25 );
        int d = (int) ( 365.25*c);
        int e = (int) ( (b - d) / 30.6001 );
        
        double dDD = b - d - (int) (30.6001*e) + f;
        int DD = (int) dDD;

        double  dhh = 24.*(dDD - DD);
        int hh = (int) dhh;
        
        double dmm = 60.*(dhh - hh);
        int mm = (int) dmm;
        
        double dss = 60.*(dmm-mm);
        int ss = (int) dss;
        //int ss = (int) (dss+0.5);
        
        int MM;
        if (e < 14) {
            MM = e - 1;
        } else {
            MM = e - 13;
        }
        
        int YYYY;
        if (MM > 2) {
            YYYY = c - 4716;
        } else {
            YYYY = c - 4715;
        }

        return new int[] {YYYY, MM, DD, hh, mm, ss};

    }
    
    public static void main(String args[]) throws Exception {

        int YYYY = Integer.parseInt(args[0]);
        int MM   = Integer.parseInt(args[1]);
        int DD   = Integer.parseInt(args[2]);
        int hh   = Integer.parseInt(args[3]);
        int mm   = Integer.parseInt(args[4]);
        int ss   = Integer.parseInt(args[5]);
        
        MJD MJD = new MJD();
        
        System.out.println("MJD for " + YYYY + "-" + MM + "-" + DD + 
                " " + hh + ":" + mm + ":" + ss + " : " +
                + MJD.greg2MJD( new int[] {YYYY, MM, DD, hh, mm, ss } ));

        int results[] = MJD.mjd2Greg(MJD.greg2MJD(new int[] {YYYY, MM, DD, hh, mm, ss } ));
        System.out.println
        ("... back to calendar : " + results[0] + " "
                + results[1] + " " + results[2] + " "
                + results[3] + " " + results[4] + " "
                + results[5]);
        
    }
    
}
