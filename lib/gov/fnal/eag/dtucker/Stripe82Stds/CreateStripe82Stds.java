package gov.fnal.eag.dtucker.Stripe82Stds;
import java.sql.*;
import java.text.*;

/**
 * This class creates a standard star catalog based upon the 
 * Stripe 82 Runs Database loads the resulting catalog into
 * a table of a SQL database.
 * @author dtucker
 *
 */
public class CreateStripe82Stds {
    
    public static void main(String[] args) throws Exception {
        
        java.util.Date d = new java.util.Date();
        System.out.println(d.toLocaleString());
        
        NumberFormat formatter1 = new DecimalFormat("0.####");
        NumberFormat formatter2 = new DecimalFormat("0.######");
        
        Class.forName("org.postgresql.Driver");
        String url = "jdbc:postgresql://sdssdp24.fnal.gov/images?user=dtucker";
        Connection db = DriverManager.getConnection(url);
        
        System.out.println("Finding min(ra), max(ra), min(dec), max(dec) for template_objects.");
        String query0 = "select min(ra), max(ra), min(dec), max(dec) from template_objects";
        Statement st0 = db.createStatement();
        ResultSet rs0 = st0.executeQuery(query0);
        rs0.next();
        double raGlobalMin  = rs0.getDouble(1);
        double raGlobalMax  = rs0.getDouble(2);
        double decGlobalMin = rs0.getDouble(3);
        double decGlobalMax = rs0.getDouble(4);
        rs0.close();
        st0.close();
        
        System.out.println("min(ra)  = " + raGlobalMin);
        System.out.println("max(ra)  = " + raGlobalMax);
        System.out.println("min(dec) = " + decGlobalMin);
        System.out.println("max(dec) = " + decGlobalMax);
        
        d = new java.util.Date();
        System.out.println(d.toLocaleString());
        
        // (status & 512) !=0         --> ok_scanline
        // objc_type = 6              --> classified as a star by photo
        // (objc_flags & 262144) = 0  --> not saturated
        // (objc_flags & 4) = 0       --> not an edge object
        // (objc_flags & 4096) = 0    --> no cosmic ray contamination
        // psf_mag_1-psf_mag_2 < 2.00 --> avoid extremely red objects
        String query1 = "SELECT ra, dec " +
        "FROM template_objects " + 
        "WHERE (ra  between ? and ?) AND " +
        "      (dec between -1.25 and  1.25) AND " +
        "      (psf_mag_2 between 15.0 and 20.0) AND " +
        "      (psf_magerr_2 < 0.03) AND " +
        "      ((status & 512) != 0) AND " +
        "      (objc_type = 6) AND " +
        "      ((objc_flags & 262144) = 0) AND " + 
        "      ((objc_flags & 4) = 0) AND " +
        "      ((objc_flags & 4096) = 0) AND " + 
        "      (psf_mag_1-psf_mag_2 < 2.00)";
        
        PreparedStatement st1 = db.prepareStatement(query1);
        
        String query2 = "SELECT psf_mag_0, psf_mag_1, psf_mag_2, psf_mag_3, psf_mag_4, " + 
        "psf_magerr_0, psf_magerr_1, psf_magerr_2, psf_magerr_3,  psf_magerr_4 " + 
        "FROM calibrated_objects " + 
        "WHERE (ra  between ? and ?) AND " +
        "      (dec between ? and ?) AND " +
        "      (psf_mag_2 between 14.5 and 20.5) AND " +
        "      (psf_magerr_2 < 0.03) AND " +
        "      ((status & 512) != 0) AND " +
        "      (objc_type = 6) AND " +
        "      ((objc_flags & 262144) = 0) AND " + 
        "      ((objc_flags & 4) = 0) AND " +
        "      ((objc_flags & 4096) = 0) AND " + 
        "      (psf_mag_1-psf_mag_2 < 2.00)";
        
        PreparedStatement st2 = db.prepareStatement(query2);
        
        int run;
        int rerun;
        double ra1;
        double ra2;
        double dec1;
        double dec2;
        double raLo;
        double raHi;
        double decLo;
        double decHi;
        double[] sum        = new double[5];
        double[] sum2       = new double[5];
        double[] psf_mag    = new double[5];
        double[] psf_magerr = new double[5];
        int[]   ntot       = new int[5]; 
        double temp1;
        double temp2;
        
        int i = 0;
        
        //double raLocalMin = raGlobalMin;
        //double raLocalMax = raLocalMin + 1.00 ;
        double raLocalMin = 25.0;
        double raLocalMax = 26.0;
        raGlobalMax = 26.0;
        
        while (raLocalMin < raGlobalMax) {
            
            st1.setDouble(1,raLocalMin);
            st1.setDouble(2,raLocalMax);
            
            ResultSet rs1 = st1.executeQuery();
            
            while (rs1.next()) {
                
                ra1  = rs1.getDouble("ra");
                dec1 = rs1.getDouble("dec");
                
                raLo  = ra1  - 0.0001;
                raHi  = ra1  + 0.0001;
                decLo = dec1 - 0.0001;
                decHi = dec1 + 0.0001;
                
                st2.setDouble(1,raLo);
                st2.setDouble(2,raHi);
                st2.setDouble(3,decLo);
                st2.setDouble(4,decHi);
                
                for (int j=0; j<5; j++) {   
                    ntot[j]       = 0; 
                    sum[j]        = 0.0;
                    sum2[j]       = 0.0;
                    psf_mag[j]    = 0.0;
                    psf_magerr[j] = 0.0;
                }
                
                ResultSet rs2 = st2.executeQuery();
                while (rs2.next()) {
                    for (int j=0; j<5; j++) {
                        String magColumnName    = "psf_mag_"+j;
                        String magerrColumnName = "psf_magerr_"+j;
                        temp1 = (double) rs2.getFloat(magColumnName);
                        temp2 = (double) rs2.getFloat(magerrColumnName);
                        if (temp2 < 0.03) {
                            sum[j]  = sum[j] + temp1;
                            sum2[j] = sum2[j] + temp1*temp1;
                            ntot[j]++;
                        }
                    }
                }                   
                rs2.close();
                
                for (int j=0; j<5; j++) {
                    psf_mag[j]    = -1000.0;
                    psf_magerr[j] = -1000.0;
                    if (ntot[j] > 0) {
                        psf_mag[j]    = sum[j]/ntot[j];
                    } 
                    if (ntot[j] > 1){
                        psf_magerr[j] = (ntot[j]*(sum2[j]/ntot[j] - psf_mag[j]*psf_mag[j])/(ntot[j]-1))/ntot[j];
                        psf_magerr[j] = Math.sqrt(psf_magerr[j]);
                    }
                }
                
                //only consider those stars with ntot(r) > 5
                // with positive psf_magerr(r), and with 
                // sigma_obs(r) < 0.05...
                if (ntot[2] > 5 && psf_magerr[2] > 0. && Math.sqrt(ntot[2])*psf_magerr[2] < 0.05 ) {
                    
                    i++;
                                        
                    String values = "'Star" + i + "', " + 
                                    formatter2.format(ra1)  + ", " + 
                                    formatter2.format(dec1) + ", " + 
                                    formatter1.format(psf_mag[0]) + ", " + 
                                    formatter1.format(psf_mag[1]) + ", " + 
                                    formatter1.format(psf_mag[2]) + ", " + 
                                    formatter1.format(psf_mag[3]) + ", " + 
                                    formatter1.format(psf_mag[4]) + ", " + 
                                    formatter1.format(psf_magerr[0]) + ", " + 
                                    formatter1.format(psf_magerr[1]) + ", " + 
                                    formatter1.format(psf_magerr[2]) + ", " + 
                                    formatter1.format(psf_magerr[3]) + ", " + 
                                    formatter1.format(psf_magerr[4]) + ", " + 
                                    ntot[0] + ", " + 
                                    ntot[1] + ", " + 
                                    ntot[2] + ", " + 
                                    ntot[3] + ", " + 
                                    ntot[4];
                    
                    if (i%1000 == 1) {
                        System.out.println(values);                
                    }
                                   
                    Statement stmt = db.createStatement();
                    stmt.executeUpdate("INSERT INTO des_stripe82_stds_test " +
                                       "VALUES (" + values + " )");
                    stmt.close();
                    
                }
                
            }
            
            rs1.close();
            
            raLocalMin = raLocalMin + 1.;
            raLocalMax = raLocalMin + 1.;
            
            d = new java.util.Date();
            System.out.println(d.toLocaleString());
            
        }
        
        
        st1.close();
        st2.close();
        
        db.close();
        
        d = new java.util.Date();
        System.out.println(d.toLocaleString());
        
        System.out.println("Finished!");
        
    }
  
}
