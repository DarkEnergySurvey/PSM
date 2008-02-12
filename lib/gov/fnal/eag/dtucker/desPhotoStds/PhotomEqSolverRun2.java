/**
 * PhotomEqSolverRun2.java   The DES Collaboration  2006
 */
package gov.fnal.eag.dtucker.desPhotoStds;

import gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolver2;
import java.sql.SQLException;
import java.util.Date;


/**
 * @author dtucker
 *
 */
public class PhotomEqSolverRun2 {


    public static void main (String[] args) {

    		System.out.println("PhotomEqSolverRun2");
    	
        System.out.print("arglist:  ");
        for (int i=0; i< args.length; i++) {
            System.out.print(args[i] + " ");
        }
        System.out.print("\n");
        
        PhotomEqSolver2 ph = new PhotomEqSolver2();

        // Process any arguments passed to the main method.
        // Assume they are ordered "url dbName username password instrument telescope nite filter ccdid magLo magHi niterations nsigma 
        //    imageType imageNameFilter psmVersion verbose".
        if (args.length > 0) {
            String url = args[0];
            ph.setUrl(url);   
            System.out.println("url="+url);
        } 
        if (args.length > 1) {
            String dbName = args[1];
            ph.setDbName(dbName);   
            System.out.println("dbName="+dbName);
        } 
        if (args.length > 2) {
            String user = args[2];
            ph.setUser(user);   
            System.out.println(user);
        } 
        if (args.length > 3) {
            String passwd = args[3];
            ph.setPasswd(passwd);   
            System.out.println(passwd);
        }
        if (args.length > 4) {
            String instrument = args[4];
            ph.setInstrument(instrument);   
            System.out.println("instrument="+instrument);
        } 
        if (args.length > 5) {
            String telescope = args[5];
            ph.setTelescope(telescope);   
            System.out.println("telescope="+telescope);
        } 
        if (args.length > 6) {
            String nite = args[6];
            ph.setNite(nite);
            System.out.println("nite="+nite);
        }
        if (args.length > 7) {
            String filter = args[7];
            ph.setFilter(filter);
            System.out.println("filter="+filter);
        } 
        if (args.length > 8) {
            int ccdid = Integer.parseInt(args[8]);
            ph.setCcdid(ccdid);
            System.out.println("ccdid="+ccdid);
        } 
        if (args.length > 9) {
            double magLo = Double.parseDouble(args[9]);
            ph.setMagLo(magLo);
            System.out.println("magLo="+magLo);
        }        
        if (args.length > 10) {
            double magHi = Double.parseDouble(args[10]);
            ph.setMagHi(magHi);
            System.out.println("magHi="+magHi);
        } 
        if (args.length > 11) {
            int niterations = Integer.parseInt(args[11]);
            ph.setNiterations(niterations);
            System.out.println("niterations="+niterations);
        } 
        if (args.length > 12) {
            double nsigma = Double.parseDouble(args[12]);
            ph.setNsigma(nsigma);   
            System.out.println("nsigma="+nsigma);
        } 
        if (args.length > 13) {
        		String imageType = args[13];
        		ph.setImageType(imageType);   
        		System.out.println("imageType="+imageType);
        } 
        if (args.length > 14) {
    			String imageNameFilter = args[14];
    			ph.setImageNameFilter(imageNameFilter);   
    			System.out.println("imageNameFilter="+imageNameFilter);
        } 
        //if (args.length > 13) {
            //String obsTable = args[13];
            //ph.setObsTable(obsTable);   
            //System.out.println("obsTable="+obsTable);
        //} 
        //if (args.length > 14) {
            //String filesTable = args[14];
            //ph.setFilesTable(filesTable);   
            //System.out.println("filesTable="+filesTable);
        //} 
        //if (args.length > 15) {
            //String stdTable = args[15];
            //ph.setStdTable(stdTable);   
            //System.out.println("stdTable="+stdTable);
        //} 
        //if (args.length > 16) {
            //String fitTable = args[16];
            //ph.setFitTable(fitTable);   
            //System.out.println("fitTable="+fitTable);
        //} 
        if (args.length > 15) {
            String runiddesc = args[15];
            ph.setRuniddesc(runiddesc);   
            System.out.println("runiddesc="+runiddesc);
        } 
        if (args.length > 16) {
            String psmVersion = args[16];
            ph.setPsmVersion(psmVersion);   
            System.out.println("psmVersion="+psmVersion);
        } 
        if (args.length > 17) {
            int verbose = Integer.parseInt(args[17]);
            ph.setVerbose(verbose);   
            System.out.println("verbose="+verbose);
        } 
       
        ph.setSqlDriver("oracle.jdbc.driver.OracleDriver");
        //ph.setUrl("jdbc:oracle:thin:@//sky.ncsa.uiuc.edu:1521/");
        //ph.setDbName("des");
        //ph.setStdTable("des_stripe82_stds_v1");
        ph.setStdTable("standard_stars");
        ph.setObsTable("OBJECTS");
        ph.setFilesTable("FILES");
        ph.setFitTable("psmfit");
        
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

    }

}
