/**
 * PhotomEqSolverTest.java The DES Collaboration  2005
 */
package gov.fnal.eag.dtucker.desPhotoStds.test;

import gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolver;
import java.sql.*;
import java.util.Date;
import junit.framework.*;
import junit.textui.TestRunner;

/**
 * @author dtucker
 *
 */
public class PhotomEqSolverTest extends TestCase {

    public static void main(String[] args) {
        TestRunner.run(suite());
    }
    
    public static Test suite() {
        return new TestSuite(PhotomEqSolverTest.class);
    }

    public void testSolveG() {
        
        // Define assert values and tolerances
        double aAssertValue = -19.900;
        double bAssertValue =   0.000;
        double kAssertValue =   0.250;
        double aAssertToler =   0.005;
        double bAssertToler =   0.005;
        double kAssertToler =   0.005;
        
        // Define database access info
        String sqlDriver  = "oracle.jdbc.driver.OracleDriver";
        //String url        = "jdbc:oracle:thin:@//talamlx.ncsa.uiuc.edu:1521/";
        String url        = "jdbc:oracle:thin:@//sky.ncsa.uiuc.edu:1521/";
        String dbName     = "des";
        String user       = "pipeline";
        String passwd     = "dc01user";

        // Define tables to be used in database
        String stdTable   = "des_stripe82_stds_v1";
        String obsTable   = "des_testobs";
        String fitTable   = "des_testfit";
        
        // Define input values for photometric solution
        double mjdLo      = 55135.0;
        double mjdHi      = 55135.4;
        String filter     = "g";
        int ccdid         = -1;
        int verbose       =  1;
        String psmVersion = "v1.0";
        

        // Fit photometric equation based upon data in stdTable and obsTable
        // and ingest results as an entry to fitTable
        PhotomEqSolver ph = new PhotomEqSolver();
       
        ph.setMjdLo(mjdLo);
        ph.setMjdHi(mjdHi);
        ph.setFilter(filter);
        ph.setCcdid(ccdid);
        ph.setVerbose(verbose);   
        ph.setPsmVersion(psmVersion);   
                
        ph.setSqlDriver(sqlDriver);
        ph.setUrl(url);
        ph.setDbName(dbName);
        ph.setUser(user);   
        ph.setPasswd(passwd);   

        ph.setStdTable(stdTable);
        ph.setObsTable(obsTable);
        ph.setFitTable(fitTable);
        
        Date date = new Date();
        ph.setDate(date);        

        try {
            ph.solve();
        } catch (ClassNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (SQLException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }        
        
        
        // Connect again to database and check whether 
        // results in the just-generated entry to testFit
        // are equal to the asserted values within the
        // tolerances define above
        try {

            Class.forName(sqlDriver);

            String fullURL = url+dbName;
            Connection db;

            db = DriverManager.getConnection(fullURL, user, passwd);
            
            Timestamp tt0 = new Timestamp(date.getTime());
            long tt0msec = tt0.getTime();
            
            String query0 = null;
            query0 = "SELECT * FROM " + fitTable;
            System.out.println("query0 = " + query0);
            
            Statement st0;
            st0 = db.createStatement();
            
            ResultSet rs0 = st0.executeQuery(query0);
            while (rs0.next()) {
                
                Timestamp tt1 = rs0.getTimestamp("timestamp");
                long tt1msec = tt1.getTime();
                long dmsec = tt1msec-tt0msec;

                // Find the entry in testFit that is within 1000 millisec
                // of the intended timestamp based upon the date variable 
                // above
                if (Math.abs(dmsec) < 1000) {

                    double a        = rs0.getDouble("a");
                    double b        = rs0.getDouble("b");
                    double k        = rs0.getDouble("k");

                    System.out.println(a + " " + b + " " + k + " " + 
                        tt1.toString() + " " + tt1.getTime() + " " + dmsec);

                    assertTrue( (a > aAssertValue-aAssertToler) && 
                        (a < aAssertValue+aAssertToler) );
                    assertTrue( (b > bAssertValue-bAssertToler) && 
                        (b < bAssertValue+bAssertToler) );
                    assertTrue( (k > kAssertValue-kAssertToler) && 
                        (k < kAssertValue+kAssertToler) );
                        
                }

            }
            
            rs0.close();
            st0.close();
            
        } catch (ClassNotFoundException e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        } catch (SQLException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }
    

    public void testSolveR() {
        
        // Define assert values and tolerances
        double aAssertValue = -20.000;
        double bAssertValue =   0.000;
        double kAssertValue =   0.150;
        double aAssertToler =   0.005;
        double bAssertToler =   0.005;
        double kAssertToler =   0.005;
        
        // Define database access info
        String sqlDriver  = "oracle.jdbc.driver.OracleDriver";
        //String url        = "jdbc:oracle:thin:@//talamlx.ncsa.uiuc.edu:1521/";
        String url        = "jdbc:oracle:thin:@//sky.ncsa.uiuc.edu:1521/";
        String dbName     = "des";
        String user       = "pipeline";
        String passwd     = "dc01user";

        // Define tables to be used in database
        String stdTable   = "des_stripe82_stds_v1";
        String obsTable   = "des_testobs";
        String fitTable   = "des_testfit";
        
        // Define input values for photometric solution
        double mjdLo      = 55135.0;
        double mjdHi      = 55135.4;
        String filter     = "r";
        int ccdid         = -1;
        int verbose       =  1;
        String psmVersion = "v1.0";
        

        // Fit photometric equation based upon data in stdTable and obsTable
        // and ingest results as an entry to fitTable
        PhotomEqSolver ph = new PhotomEqSolver();
       
        ph.setMjdLo(mjdLo);
        ph.setMjdHi(mjdHi);
        ph.setFilter(filter);
        ph.setCcdid(ccdid);
        ph.setVerbose(verbose);   
        ph.setPsmVersion(psmVersion);   
                
        ph.setSqlDriver(sqlDriver);
        ph.setUrl(url);
        ph.setDbName(dbName);
        ph.setUser(user);   
        ph.setPasswd(passwd);   

        ph.setStdTable(stdTable);
        ph.setObsTable(obsTable);
        ph.setFitTable(fitTable);
        
        Date date = new Date();
        ph.setDate(date);        

        try {
            ph.solve();
        } catch (ClassNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (SQLException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }        
        
        
        // Connect again to database and check whether 
        // results in the just-generated entry to testFit
        // are equal to the asserted values within the
        // tolerances define above
        try {

            Class.forName(sqlDriver);

            String fullURL = url+dbName;
            Connection db;

            db = DriverManager.getConnection(fullURL, user, passwd);
            
            Timestamp tt0 = new Timestamp(date.getTime());
            long tt0msec = tt0.getTime();
            
            String query0 = null;
            query0 = "SELECT * FROM " + fitTable;
            System.out.println("query0 = " + query0);
            
            Statement st0;
            st0 = db.createStatement();
            
            ResultSet rs0 = st0.executeQuery(query0);
            while (rs0.next()) {
                
                Timestamp tt1 = rs0.getTimestamp("timestamp");
                long tt1msec = tt1.getTime();
                long dmsec = tt1msec-tt0msec;

                // Find the entry in testFit that is within 1000 millisec
                // of the intended timestamp based upon the date variable 
                // above
                if (Math.abs(dmsec) < 1000) {

                    double a        = rs0.getDouble("a");
                    double b        = rs0.getDouble("b");
                    double k        = rs0.getDouble("k");

                    System.out.println(a + " " + b + " " + k + " " + 
                        tt1.toString() + " " + tt1.getTime() + " " + dmsec);

                    assertTrue( (a > aAssertValue-aAssertToler) && 
                        (a < aAssertValue+aAssertToler) );
                    assertTrue( (b > bAssertValue-bAssertToler) && 
                        (b < bAssertValue+bAssertToler) );
                    assertTrue( (k > kAssertValue-kAssertToler) && 
                        (k < kAssertValue+kAssertToler) );
                        
                }

            }
            
            rs0.close();
            st0.close();
            
        } catch (ClassNotFoundException e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        } catch (SQLException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    
    }
    

    public void testSolveI() {
        
        // Define assert values and tolerances
        double aAssertValue = -20.100;
        double bAssertValue =   0.000;
        double kAssertValue =   0.100;
        double aAssertToler =   0.005;
        double bAssertToler =   0.005;
        double kAssertToler =   0.005;
        
        // Define database access info
        String sqlDriver  = "oracle.jdbc.driver.OracleDriver";
        //String url        = "jdbc:oracle:thin:@//talamlx.ncsa.uiuc.edu:1521/";
        String url        = "jdbc:oracle:thin:@//sky.ncsa.uiuc.edu:1521/";
        String dbName     = "des";
        String user       = "pipeline";
        String passwd     = "dc01user";

        // Define tables to be used in database
        String stdTable   = "des_stripe82_stds_v1";
        String obsTable   = "des_testobs";
        String fitTable   = "des_testfit";
        
        // Define input values for photometric solution
        double mjdLo      = 55135.0;
        double mjdHi      = 55135.4;
        String filter     = "i";
        int ccdid         = -1;
        int verbose       =  1;
        String psmVersion = "v1.0";
        

        // Fit photometric equation based upon data in stdTable and obsTable
        // and ingest results as an entry to fitTable
        PhotomEqSolver ph = new PhotomEqSolver();
       
        ph.setMjdLo(mjdLo);
        ph.setMjdHi(mjdHi);
        ph.setFilter(filter);
        ph.setCcdid(ccdid);
        ph.setVerbose(verbose);   
        ph.setPsmVersion(psmVersion);   
                
        ph.setSqlDriver(sqlDriver);
        ph.setUrl(url);
        ph.setDbName(dbName);
        ph.setUser(user);   
        ph.setPasswd(passwd);   

        ph.setStdTable(stdTable);
        ph.setObsTable(obsTable);
        ph.setFitTable(fitTable);
        
        Date date = new Date();
        ph.setDate(date);        

        try {
            ph.solve();
        } catch (ClassNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (SQLException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }        
        
        
        // Connect again to database and check whether 
        // results in the just-generated entry to testFit
        // are equal to the asserted values within the
        // tolerances define above
        try {

            Class.forName(sqlDriver);

            String fullURL = url+dbName;
            Connection db;

            db = DriverManager.getConnection(fullURL, user, passwd);
            
            Timestamp tt0 = new Timestamp(date.getTime());
            long tt0msec = tt0.getTime();
            
            String query0 = null;
            query0 = "SELECT * FROM " + fitTable;
            System.out.println("query0 = " + query0);
            
            Statement st0;
            st0 = db.createStatement();
            
            ResultSet rs0 = st0.executeQuery(query0);
            while (rs0.next()) {
                
                Timestamp tt1 = rs0.getTimestamp("timestamp");
                long tt1msec = tt1.getTime();
                long dmsec = tt1msec-tt0msec;

                // Find the entry in testFit that is within 1000 millisec
                // of the intended timestamp based upon the date variable 
                // above
                if (Math.abs(dmsec) < 1000) {

                    double a        = rs0.getDouble("a");
                    double b        = rs0.getDouble("b");
                    double k        = rs0.getDouble("k");

                    System.out.println(a + " " + b + " " + k + " " + 
                        tt1.toString() + " " + tt1.getTime() + " " + dmsec);

                    assertTrue( (a > aAssertValue-aAssertToler) && 
                        (a < aAssertValue+aAssertToler) );
                    assertTrue( (b > bAssertValue-bAssertToler) && 
                        (b < bAssertValue+bAssertToler) );
                    assertTrue( (k > kAssertValue-kAssertToler) && 
                        (k < kAssertValue+kAssertToler) );
                        
                }

            }
            
            rs0.close();
            st0.close();
            
        } catch (ClassNotFoundException e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        } catch (SQLException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }
    

    public void testSolveZ() {
        
        // Define assert values and tolerances
        double aAssertValue = -20.200;
        double bAssertValue =   0.000;
        double kAssertValue =   0.050;
        double aAssertToler =   0.005;
        double bAssertToler =   0.005;
        double kAssertToler =   0.005;
        
        // Define database access info
        String sqlDriver  = "oracle.jdbc.driver.OracleDriver";
        //String url        = "jdbc:oracle:thin:@//talamlx.ncsa.uiuc.edu:1521/";
        String url        = "jdbc:oracle:thin:@//sky.ncsa.uiuc.edu:1521/";
        String dbName     = "des";
        String user       = "pipeline";
        String passwd     = "dc01user";

        // Define tables to be used in database
        String stdTable   = "des_stripe82_stds_v1";
        String obsTable   = "des_testobs";
        String fitTable   = "des_testfit";
        
        // Define input values for photometric solution
        double mjdLo      = 55135.0;
        double mjdHi      = 55135.4;
        String filter     = "z";
        int ccdid         = -1;
        int verbose       =  1;
        String psmVersion = "v1.0";
        

        // Fit photometric equation based upon data in stdTable and obsTable
        // and ingest results as an entry to fitTable
        PhotomEqSolver ph = new PhotomEqSolver();
       
        ph.setMjdLo(mjdLo);
        ph.setMjdHi(mjdHi);
        ph.setFilter(filter);
        ph.setCcdid(ccdid);
        ph.setVerbose(verbose);   
        ph.setPsmVersion(psmVersion);   
                
        ph.setSqlDriver(sqlDriver);
        ph.setUrl(url);
        ph.setDbName(dbName);
        ph.setUser(user);   
        ph.setPasswd(passwd);   

        ph.setStdTable(stdTable);
        ph.setObsTable(obsTable);
        ph.setFitTable(fitTable);
        
        Date date = new Date();
        ph.setDate(date);        

        try {
            ph.solve();
        } catch (ClassNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (SQLException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }        
        
        
        // Connect again to database and check whether 
        // results in the just-generated entry to testFit
        // are equal to the asserted values within the
        // tolerances define above
        try {

            Class.forName(sqlDriver);

            String fullURL = url+dbName;
            Connection db;

            db = DriverManager.getConnection(fullURL, user, passwd);
            
            Timestamp tt0 = new Timestamp(date.getTime());
            long tt0msec = tt0.getTime();
            
            String query0 = null;
            query0 = "SELECT * FROM " + fitTable;
            System.out.println("query0 = " + query0);
            
            Statement st0;
            st0 = db.createStatement();
            
            ResultSet rs0 = st0.executeQuery(query0);
            while (rs0.next()) {
                
                Timestamp tt1 = rs0.getTimestamp("timestamp");
                long tt1msec = tt1.getTime();
                long dmsec = tt1msec-tt0msec;

                // Find the entry in testFit that is within 1000 millisec
                // of the intended timestamp based upon the date variable 
                // above
                if (Math.abs(dmsec) < 1000) {

                    double a        = rs0.getDouble("a");
                    double b        = rs0.getDouble("b");
                    double k        = rs0.getDouble("k");

                    System.out.println(a + " " + b + " " + k + " " + 
                        tt1.toString() + " " + tt1.getTime() + " " + dmsec);

                    assertTrue( (a > aAssertValue-aAssertToler) && 
                        (a < aAssertValue+aAssertToler) );
                    assertTrue( (b > bAssertValue-bAssertToler) && 
                        (b < bAssertValue+bAssertToler) );
                    assertTrue( (k > kAssertValue-kAssertToler) && 
                        (k < kAssertValue+kAssertToler) );
                        
                }

            }
            
            rs0.close();
            st0.close();
            
        } catch (ClassNotFoundException e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        } catch (SQLException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }
    

}
