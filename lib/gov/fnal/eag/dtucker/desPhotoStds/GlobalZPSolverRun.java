/**
 * PhotomEqSolverRun2.java   The DES Collaboration  2006
 */
package gov.fnal.eag.dtucker.desPhotoStds;

import gov.fnal.eag.dtucker.desPhotoStds.GlobalZPSolver;
import java.sql.SQLException;
import java.util.Date;


/**
 * @author dtucker
 *
 */
public class GlobalZPSolverRun {


    public static void main (String[] args) {

    		System.out.println("GlobalZPSolverRun");
    	
        System.out.print("arglist:  ");
        for (int i=0; i< args.length; i++) {
            System.out.print(args[i] + " ");
        }
        System.out.print("\n");
        
        GlobalZPSolver globalZP = new GlobalZPSolver();

        if (args.length > 0) {
            String starMatchFileName = args[0];
            globalZP.setStarMatchFileName(starMatchFileName);   
            System.out.println("starMatchFileName="+starMatchFileName);
        } 
        if (args.length > 1) {
            String outputFileName = args[1];
            globalZP.setOutputFileName(outputFileName);   
            System.out.println("outputFileName="+outputFileName);
        } 
         if (args.length > 2) {
            int verbose = Integer.parseInt(args[2]);
            globalZP.setVerbose(verbose);   
            System.out.println("verbose="+verbose);
        } 

        Date date = new Date();
        globalZP.setDate(date);
        
        try {
            globalZP.solve();
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
