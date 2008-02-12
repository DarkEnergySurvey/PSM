/**
 * PhotomEqSolverRun2.java   The DES Collaboration  2006
 */
package gov.fnal.eag.dtucker.desPhotoStds;

import gov.fnal.eag.dtucker.desPhotoStds.StarFlatSolver;
import java.sql.SQLException;
import java.util.Date;


/**
 * @author dtucker
 *
 */
public class StarFlatSolverRun {


    public static void main (String[] args) {

    		System.out.println("StarFlatSolverRun");
    	
        System.out.print("arglist:  ");
        for (int i=0; i< args.length; i++) {
            System.out.print(args[i] + " ");
        }
        System.out.print("\n");
        
        StarFlatSolver starFlat = new StarFlatSolver();

        // Process any arguments passed to the main method.
        // Assume they are ordered "url dbName username password instrument telescope nite filter ccdid magLo magHi niterations nsigma 
        //    imageType imageNameFilter psmVersion verbose".
        if (args.length > 0) {
            String url = args[0];
            starFlat.setUrl(url);   
            System.out.println("url="+url);
        } 
        if (args.length > 1) {
            String dbName = args[1];
            starFlat.setDbName(dbName);   
            System.out.println("dbName="+dbName);
        } 
        if (args.length > 2) {
            String user = args[2];
            starFlat.setUser(user);   
            System.out.println(user);
        } 
        if (args.length > 3) {
            String passwd = args[3];
            starFlat.setPasswd(passwd);   
            System.out.println(passwd);
        }
        if (args.length > 4) {
            String instrument = args[4];
            starFlat.setInstrument(instrument);   
            System.out.println("instrument="+instrument);
        } 
        if (args.length > 5) {
            String telescope = args[5];
            starFlat.setTelescope(telescope);   
            System.out.println("telescope="+telescope);
        } 
        if (args.length > 6) {
            String nite = args[6];
            starFlat.setNite(nite);
            System.out.println("nite="+nite);
        }
        if (args.length > 7) {
            String filter = args[7];
            starFlat.setFilter(filter);
            System.out.println("filter="+filter);
        } 
        if (args.length > 8) {
            int ccdid = Integer.parseInt(args[8]);
            starFlat.setCcdid(ccdid);
            System.out.println("ccdid="+ccdid);
        } 
        if (args.length > 9) {
            double magLo = Double.parseDouble(args[9]);
            starFlat.setMagLo(magLo);
            System.out.println("magLo="+magLo);
        }        
        if (args.length > 10) {
            double magHi = Double.parseDouble(args[10]);
            starFlat.setMagHi(magHi);
            System.out.println("magHi="+magHi);
        } 
        if (args.length > 11) {
            int niterations = Integer.parseInt(args[11]);
            starFlat.setNiterations(niterations);
            System.out.println("niterations="+niterations);
        } 
        if (args.length > 12) {
            double nsigma = Double.parseDouble(args[12]);
            starFlat.setNsigma(nsigma);   
            System.out.println("nsigma="+nsigma);
        } 
        if (args.length > 13) {
        		String imageType = args[13];
        		starFlat.setImageType(imageType);   
        		System.out.println("imageType="+imageType);
        } 
        if (args.length > 14) {
    			String imageNameFilter = args[14];
    			starFlat.setImageNameFilter(imageNameFilter);   
    			System.out.println("imageNameFilter="+imageNameFilter);
        } 
        //if (args.length > 13) {
            //String obsTable = args[13];
            //starFlat.setObsTable(obsTable);   
            //System.out.println("obsTable="+obsTable);
        //} 
        //if (args.length > 14) {
            //String filesTable = args[14];
            //starFlat.setFilesTable(filesTable);   
            //System.out.println("filesTable="+filesTable);
        //} 
        //if (args.length > 15) {
            //String stdTable = args[15];
            //starFlat.setStdTable(stdTable);   
            //System.out.println("stdTable="+stdTable);
        //} 
        //if (args.length > 16) {
            //String fitTable = args[16];
            //starFlat.setFitTable(fitTable);   
            //System.out.println("fitTable="+fitTable);
        //} 
        if (args.length > 15) {
            String psmVersion = args[15];
            starFlat.setPsmVersion(psmVersion);   
            System.out.println("psmVersion="+psmVersion);
        } 
        if (args.length > 16) {
            int verbose = Integer.parseInt(args[16]);
            starFlat.setVerbose(verbose);   
            System.out.println("verbose="+verbose);
        } 
       
        starFlat.setSqlDriver("oracle.jdbc.driver.OracleDriver");
        //starFlat.setUrl("jdbc:oracle:thin:@//sky.ncsa.uiuc.edu:1521/");
        //starFlat.setDbName("des");
        //starFlat.setStdTable("des_stripe82_stds_v1");
        starFlat.setStdTable("standard_stars");
        starFlat.setObsTable("OBJECTS");
        starFlat.setFilesTable("FILES");
        starFlat.setFitTable("psmfit");
        
        Date date = new Date();
        starFlat.setDate(date);
        
        try {
            starFlat.solve();
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

    }

}
