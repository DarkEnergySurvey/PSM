package gov.fnal.eag.dtucker.desPhotoStds;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Random;

public class TestObsSQLTable {
    
    private float[] a    = new float[5];
    private float[] k    = new float[5];
    private float[] rms  = new float[5];
    
    public static void main(String[] args) throws ClassNotFoundException {
        
        TestObsSQLTable tt = new TestObsSQLTable();
        Class.forName("org.postgresql.Driver");
        String url = "jdbc:postgresql://sdssdp24.fnal.gov/images?user=dtucker";
        
        Connection db;
        try {
            db = DriverManager.getConnection(url);
            tt.createTestObsSQLTable(db);
            tt.populateTestObsSQLTable(db);
            db.close();
        } catch (SQLException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }        
        
    }
    
    public void createTestObsSQLTable (Connection db) throws ClassNotFoundException, SQLException {
        
        String command = "CREATE TABLE des_testobs ( " +
        "ra  double precision, "        +
        "dec double precision, "        +
        "mjd double precision, "        +
        "ccd int, "           +
        "exptime real, "      + 
        "airmass real, "      +
        "filter varchar(2), " +
        "counts1 real, "      +
        "counts2 real, "      +
        "counts3 real, "      +
        "saturated int, "     + 
        "a float, "           +
        "k float, "           +
        "rms float, "         +
        "stdmag float, "      +
        "instmag float );";
        
        Statement stmt = db.createStatement();
        stmt.executeUpdate(command);
        stmt.close();
        
        System.out.println("Finished with createTestObsSQLTable");

    }
    
    public void populateTestObsSQLTable (Connection db) throws SQLException {

        //mjd0 = Oct 31, 2009
        double mjd0 = 55135.0;
        
        float[] stdmag   = new float[5];
        float[] instmag  = new float[5];
        float[] counts   = new float[5];
        float[] satLimit = new float[5];
        int[] saturated  = new int[5];
        String[] filter = new String[5];
        
        a[0] = -19.8f;
        a[1] = -19.9f;
        a[2] = -20.0f;
        a[3] = -20.1f;
        a[4] = -20.2f;
        
        k[0] = 0.60f;
        k[1] = 0.25f;
        k[2] = 0.15f;
        k[3] = 0.10f;
        k[4] = 0.05f;
        
        rms[0] = 0.040f;
        rms[1] = 0.020f;
        rms[2] = 0.020f;
        rms[3] = 0.020f;
        rms[4] = 0.020f;
        
        filter[0] = "u";
        filter[1] = "g";
        filter[2] = "r";
        filter[3] = "i";
        filter[4] = "z";
        
        int ccdid = 29;
        float exptime = 100.0f;
        for (int i=0; i<5; i++) {
            satLimit[i] = (float) (exptime*Math.pow(10.0,-0.4*(15.5+a[i])));
        }

        int ientry = 0;
        
        Random generator = new Random();
        
        String query1 = "SELECT * FROM des_stripe82_stds " +
        "WHERE (radeg between 15.0 and 16.0) AND " +
        "      (decdeg between -0.50 and  0.50) AND " +
        "      (stdmag_r between 15.0 and 20.0)";
        
        Statement st1 = db.createStatement();
        ResultSet rs1 = st1.executeQuery(query1);
        
        while (rs1.next()) {
            
            double ra  = rs1.getDouble("radeg");
            double dec = rs1.getDouble("decdeg");
            stdmag[0] = rs1.getFloat("stdmag_u");
            stdmag[1] = rs1.getFloat("stdmag_g");
            stdmag[2] = rs1.getFloat("stdmag_r");
            stdmag[3] = rs1.getFloat("stdmag_i");
            stdmag[4] = rs1.getFloat("stdmag_z");
            
            for (int j=0; j<5; j++) {
                
                float X = 1.0f + j*0.2f;
                double mjd = mjd0 + (j+3.0)/24.0;
                
                for (int i=0; i<5; i++) {
                    saturated[i] = 0;
                    if (stdmag[i] > 0.) {
                        instmag[i] = (float) ((stdmag[i]+a[i]+k[i]*X)+rms[i]*generator.nextGaussian());
                        counts[i] = (float) (exptime * Math.pow(10.0,-0.4*instmag[i]));
                        if (counts[i] > satLimit[i]) {
                            saturated[i] = 1;
                        }
                    } else {
                        instmag[i] = -1000.0f;
                        counts[i]  = -1000.0f;
                    }

                    String values = ra  + ", " + dec + ", " + mjd + ", " + 
                                    ccdid + ", " + exptime + ", " + X + ", '" + filter[i] + "', " +
                                    counts[i] + ", " + counts[i] + ", " + counts[i] + ", " + 
                                    saturated[i] + ", " + 
                                    a[i] + ", " + k[i] + ", " + rms[i] + ", " + 
                                    stdmag[i] + ", " + instmag[i];
                                       
                    Statement stmt = db.createStatement();
                    stmt.executeUpdate("INSERT INTO des_testobs " + 
                            "VALUES (" + values + " )");
                    stmt.close();
                    
                    
                    if (ientry%1000 == 1) {
                        System.out.println(values);                
                    }
                    ientry++;
                    
                }
                
            }      
            
        }
        
        rs1.close();
        st1.close();
       
        System.out.println("Finished with populateTestObsSQLTable");
        
    }
    
}
















