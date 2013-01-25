package gov.fnal.eag.dtucker.desPhotoStds;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.*;
import java.util.ArrayList;
import java.util.Date;
import java.util.StringTokenizer;

import nom.tam.fits.*;
import nom.tam.util.BufferedFile;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYBoxAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBubbleRenderer;
import org.jfree.chart.renderer.xy.XYDotRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.TextAnchor;


import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

// import gov.fnal.eag.dtucker.util.MJD;

/**
 * This class defines methods for solving the photometric zeropoints 
 * (a_1, a_2, ..., a_N), the instrumental color term coefficients
 * (b_1, b_2, ..., b_N), and the first-order extinction (k) for a 
 * given filter for a given night by fitting the following equation 
 * for standard star observations: 
 *    m_inst-m_std = a_1 + ... a_N + 
 *                   b_1*(stdColor-stdColor0) + ... + 
 *                   b_N*(stdColor-stdColor0) + kX,
 * where m_inst is the instrumental (observed) mag of a standard 
 * star, m_std is the known calibrated mag of the standard star, 
 * stdColor is the known calibrated color of the standard star
 * (e.g., its g-r color), stdColor0 is a zeropoint constant for  
 * the standard color and X is the airmass of the observation.
 * 
 * For a camera with a single CCD, the above equation reduces to 
 * the following, simpler form:
 *    m_inst-m_std = a + b*(stdColor-stdColor0) + kX
 *
 * For the explicit case of the g filter and a g-r color,
 * this single-CCD example looks like this:
 *   g_inst-g_std = a + b*( (g-r) - (g-r)_0 ) + kX
 * 
 * @author dtucker
 * 
 */
public class PhotomEqSolverYear1 {

	// Instance variables dealing with the SQL database
	private String sqlDriver = "oracle.jdbc.driver.OracleDriver";
	private String url = "jdbc:oracle:thin:@charon.ncsa.uiuc.edu:1521:";
	private String dbName = "des";
	private String user = "dummy";
	private String passwd = "dummy";
	private String obsTable = "OBJECTS";
	private String imageTable = "IMAGE";
	private String stdTable = "standard_stars_all";
	private int standard_set_in = 3;
	private String fitTable = "psmfit";

	// Instance variables dealing with this execution of the Photometric 
	// Standards Module
	private String[] filterList = { "u", "g", "r", "i", "z", "Y" };
	private int verbose = 0;
	private String psmVersion = "v_Year1";
	private Date date = new Date();
	private boolean skipQAPlots = false;

	// Instance variables dealing with the observed data to be 
	// calibrated
	private String project = "DES";
	private double mjdLo =  999999; 
	private double mjdHi = -999999; 
	private int ccdid = 0;     // ccdid=0 means all ccds
	private String filter = "r";
	private String nite = "2005oct09";
	private double magLo = 15.0;
	private double magHi = 18.0;
	private String imageType = "remap";
	private String imageNameFilter = "%";
	private String imageidExcludeList = "";
	private String run = "%";
	private String magType = "mag_psf";

	// Instance variables dealing with the fit
	private int niterations = 3; // for outlier removal in fit
	private double nsigma = 2.5; // for outlier removal in fit
	private boolean updateDB = false; //set to true if updating fitTable directly
	private boolean useOnlyCurrentObjects = false; //set to true if only using objects from the OBJECTS_CURRENT table
	
	private double stdColor0   = 0.00;   //reference or fiducial standard color
	private double stdColorLo  = -1000.; //lower end of range of stdColors to consider in solution
	private double stdColorHi  =  1000.; //upper end of range of stdColors to consider in solution
	private double baseMagErr = 0.02; // minimum error assoc'd with the
                                      // observed mag
	
	private boolean asolve      = true; // we don't currently use
	private double  adefault    = 0.00; // these a term private 
	private double  adefaultErr = 0.01; // variables except for PSMFIT

	private boolean bsolve      = false;
	private double  bdefault    = 0.00;
	private double  bdefaultErr = 0.01;
	private ColorTermCoeffs colorTermCoeffs = new ColorTermCoeffs();
	
	private boolean ksolve      = false;
	private double  kdefault    = 0.18;
	private double  kdefaultErr = 0.02;
	
	
	public void solve() throws Exception {

		if (verbose > 0) {
			System.out.println("\n\nPhotomEqSolverYear1");
		}
			
		if (verbose > 0) {
			System.out.println("");
			System.out.println("The beginning...");
			System.out.println("");
		}

		int filterIndex = -1;
		for (int j = 0; j < filterList.length; j++) {
			if (filter.equals(filterList[j])) {
				filterIndex = j;
				break;
			}
		}
		if (verbose > 1) {
			System.out.println("filterIndex for filter " + filter + " is "
					+ filterIndex + ".");
			System.out.println("");
		}
		if (filterIndex < 0) {
			System.out.println("STATUS5BEG ** Incompatible filter index.  Throwing Exception! ** STATUS5END");
			System.out.println("");
			throw new Exception();
		}

		//Parameters needed for determination of b terms...
		String cFilter = "";
		String stdColorName = "";
		if (filter.equals("u")) {
			cFilter = "g";
			stdColorName = "u-g";
		} else if (filter.equals("g")) {
			cFilter = "r";
			stdColorName = "g-r";
		} else if (filter.equals("r")) {
			cFilter = "g";
			stdColorName = "g-r";
		} else if (filter.equals("i")) {
			cFilter = "z";
			stdColorName = "i-z";
		} else if (filter.equals("z")) {
			cFilter = "i";
			stdColorName = "i-z";
		} else if (filter.equals("Y")) {
			cFilter = "z";
			stdColorName = "z-Y";
		}

		if (verbose > 1) {
			System.out.println("filter is " + filter + " cFilter is " + cFilter);
		}
		int cFilterIndex = -1;
		for (int j = 0; j < filterList.length; j++) {
			if (cFilter.equals(filterList[j])) {
				cFilterIndex = j;
				break;
			}
		}
		if (verbose > 1) {
			System.out.println("cFilterIndex for cFilter " + cFilter + " is "
					+ cFilterIndex + ".");
			System.out.println("");
		}
		if (cFilterIndex < 0) {
			System.out.println("STATUS5BEG ** Incompatible cFilter index.  Throwing Exception! ** STATUS5END");
			System.out.println("");
			throw new Exception();
		}
		
		// Establish connection to database
		if (verbose > 0) {
			System.out.println("Establishing connection to database.");
			System.out.println("");
		}
		Class.forName(sqlDriver);
		String fullURL = url + dbName;
		Connection db = DriverManager.getConnection(fullURL, user, passwd);
		
		String query0;
		if (useOnlyCurrentObjects) {
			query0 = "SELECT * FROM table(fPhotoStdsMatch2(" + "'" + imageType
			    + "', '" + imageNameFilter + "', " + "'" + nite + "', '"
				+ filter + "', " + ccdid + ", " + magLo + ", " + magHi + ", '"
				+ run + "', '" + project + "', 'CURRENT', '" + stdTable + "'," + standard_set_in + "))";
		} else {
			query0 = "SELECT * FROM table(fPhotoStdsMatch2(" + "'" + imageType
		    + "', '" + imageNameFilter + "', " + "'" + nite + "', '"
			+ filter + "', " + ccdid + ", " + magLo + ", " + magHi + ", '"
			+ run + "', '" + project + "', 'ALL', '" + stdTable + "'," + standard_set_in + "))";
		}
		if (verbose > 1) {
			System.out.println("query0 = " + query0);
			System.out.println("");
		}

		//Update due to new standard star table and fPhotoStdsMatch2:
		//String query1 = "SELECT * FROM " + stdTable
		//+ " WHERE standard_star_id = ?";
		String query1 = "SELECT * FROM " + stdTable
		+ " WHERE id = ?";
		PreparedStatement st1 = db.prepareStatement(query1);
		if (verbose > 1) {
			System.out.println("query1 = " + query1);
			System.out.println("");
		}

		// Might be interesting to grab all the mag apertures and the mag_psf to
		//  perform a curve-of-growth analysis/QA plot; perhaps in a later version
		//  of this code...
		String magErrType = "MAGERR" + magType.substring(3);
		//String query2 = "SELECT i.CCD, o." + magType + ", o." + magErrType + ", i.id, i.EXPOSUREID, o.X_IMAGE, o.Y_IMAGE FROM " + imageTable + " i, "
		//		+ obsTable + " o WHERE o.imageid=i.id AND o.object_id = ?";
		String query2 = "SELECT i.CCD, o." + magType + ", o." + magErrType + ", i.id, i.EXPOSUREID, o.XWIN_IMAGE, o.YWIN_IMAGE FROM " + imageTable + " i, "
		+ obsTable + " o WHERE o.imageid=i.id AND o.object_id = ?";
		PreparedStatement st2 = db.prepareStatement(query2);
		if (verbose > 1) {
			System.out.println("query2 = " + query2);
			System.out.println("");
		}

		String query3 = "SELECT MJD_OBS FROM EXPOSURE WHERE id = ?";
		PreparedStatement st3 = db.prepareStatement(query3);
		if (verbose > 1) {
			System.out.println("query3 = " + query3);
			System.out.println("");
		}

		String query4 = "SELECT max(psmfit_id) FROM " + fitTable;
		if (verbose > 1) {
			System.out.println("query4 = " + query4);
			System.out.println("");
		}
		
		
		// Create array list of image id's to be excluded from the fit...
		ArrayList imageidExcludeArrayList = new ArrayList();
		StringTokenizer st = new StringTokenizer(imageidExcludeList,",");
		int nTokens = st.countTokens();
		System.out.println("nTokens=" + nTokens);
		for (int i=0; i<nTokens; i++) {
			//int imageid2Exclude = Integer.parseInt((st.nextToken()).trim());
			long imageid2Exclude = Long.parseLong((st.nextToken()).trim());
			if (verbose > 1) {	
				System.out.println(i + "\t" + imageid2Exclude + " added to imageid exclude list");
			}
			if (imageidExcludeArrayList.contains(imageid2Exclude) == false) {
				//imageidExcludeArrayList.add(new Integer(imageid2Exclude));
				imageidExcludeArrayList.add(new Long(imageid2Exclude));
			}
		}		
		
		
		// Prepare to loop through all observations of standard stars obtained that night...
		
		if (verbose > 1) {
			System.out
					.println("Point  ccd_number  image_id  exposure_id  mjd_obs  object_id   x_image   y_image   standard_star_id  stdmag["
							+ filter + "]  stdmag[" + cFilter + "]  instmag  instmagErr  airmass  fieldName");
		}
		int i = 0;
		double[] stdmag = new double[filterList.length];
		double[] stdmagerr = new double[filterList.length];

		// Set up an array of ccd IDs...
		int ccdIdArray[] = new int[100];
		// iccd is the index of the ccdIdArray, nccd is the total 
		// number of ccds in the fit...
		int iccd = 0;
		int nccd = 0;
		// Initialize ccdIdArray...
		for (iccd=0; iccd<100; iccd++) {
			ccdIdArray[iccd] = -1;
		}

		// These variables are used later for the mag residuals vs focal plane position QA plot.
		// xCenter and yCenter denotes the position of the center of the focal plane in pixel coordinates
		// relative to the pixel coordinate system internal to a CCD.
		// nAxis1 and nAxis2 are the dimensions of the CCD in pixel coordinates.
		// (x,y)Axis(Min,Max) denote the min, max coordinates in pixels for the entire focal plane. 
		double[] xCenter  = new double[100];
		double[] yCenter   = new double[100];
		int nAxis1[] = new int[100];
		int nAxis2[] = new int[100];
 		double xAxisMin = 0.;
		double xAxisMax = 0.;
		double yAxisMin = 0.;
		double yAxisMax = 0.;
		
		Statement st0 = db.createStatement();
		ResultSet rs0 = st0.executeQuery(query0);
		ArrayList[] mStdStarList = new ArrayList[100];

		while (rs0.next()) {

			long object_id = rs0.getLong("object_id");
			int standard_star_id = rs0.getInt("standard_star_id");

			st1.setInt(1, standard_star_id);
			ResultSet rs1 = st1.executeQuery();
			rs1.next();

			String stdStarName = (String) rs1.getString("name");
			String fieldName = (String) rs1.getString("fieldname");
			for (int j = 0; j < filterList.length; j++) {
				String magColumnName = "stdmag_" + filterList[j];
				String magerrColumnName = "stdmagerr_" + filterList[j];
				stdmag[j] = (double) rs1.getFloat(magColumnName);
				stdmagerr[j] = (double) rs1.getFloat(magerrColumnName);
			}

			rs1.close();

			// If the stdmag or the stdmagerr for the filter being 
			// solved for doesn't make sense, then skip this std 
			// star...
			if (stdmag[filterIndex] < 0. || stdmag[filterIndex] > 30. || stdmagerr[filterIndex] < 0.) {
				continue;
			}
			// If the stdmag or the stdmagerr for the cFilter for 
			// the filter being solved for doesn't make sense, then 
			// skip this std star...
			if (stdmag[cFilterIndex] < 0. || stdmag[cFilterIndex] > 30. || stdmagerr[cFilterIndex] < 0.) {
				continue;
			}
			
			// If the stdmagerr in the filter being solved for is 
			// greater than 0.015mag, then skip this std star...
			if (stdmagerr[filterIndex] > 0.015) {
				continue;
			}
			
			//double totMagErr = Math.sqrt(baseMagErr*baseMagErr + stdmagerr[filterIndex]*stdmagerr[filterIndex]);
			double totMagErr = baseMagErr;
			
			double exptime = (double) rs0.getFloat("exptime");
			double zeropoint = (double) rs0.getFloat("zeropoint");
			double airmass = (double) rs0.getFloat("airmass");


			// Find on which CCD, image, and exposure this star lies, and find this star's 
			//  instrumental mag, instrumental mag err, and (x,y)-position on the CCD...
			st2.setLong(1, object_id);
			ResultSet rs2 = st2.executeQuery();
			rs2.next();
			int ccd_number = (int) rs2.getInt(1);
			double instmag = (double) rs2.getFloat(2);
			double instmagErr = (double) rs2.getFloat(3);
			long image_id = (long) rs2.getInt(4);
			long exposure_id = (long) rs2.getInt(5);
			double x_image = (double) rs2.getDouble(6);
			double y_image = (double) rs2.getDouble(7);
			rs2.close();
			
			//if (ccd_number == 60) {
			//	continue;
			//}
			
			
			// If this star's instrumental mag is too small or too large, skip it
			if (instmag <= 0.0 || instmag >= 99.) {
				continue;
			}
			// If the error in this star's instrumental mag is too small or too large, skip it
			if (instmagErr <= 0.0 || instmagErr >= 0.20) {
				continue;
			}
			
			
			// If this star lies in one of the excluded imageid's, skip it
			if (imageidExcludeArrayList.contains(image_id) == true) {
				if (verbose > 0) {
					System.out.println("imageid " + image_id + " is part of imageid exclude list...  skipping... ");
				}
				continue;
			}

			// If this star's standard color falls outside the standard color range, skip it..
			double stdColor = 0.;
			if (filter.equals("u")) {
				stdColor = stdmag[0]-stdmag[1];
			} else if (filter.equals("g") || filter.equals("r")) {
				stdColor = stdmag[1]-stdmag[2];
			} else if (filter.equals("i") || filter.equals("z")) {
				stdColor = stdmag[3]-stdmag[4];
			} else if (filter.equals("Y")) {
				stdColor = stdmag[4]-stdmag[5];
			}
			if (stdColor < stdColorLo || stdColor > stdColorHi) {
				continue;
			}
			
			// Calculate the deltamag of this observation of this star
			instmag = instmag - zeropoint;
			if (exptime > 0.) {
				instmag = instmag + 2.5 * 0.4342944819 * Math.log(exptime);
			}
			double deltamag = instmag - stdmag[filterIndex];			

			// If available, find the mjd_obs of the exposure from which this observation came:
			double mjd_obs = 0.;
			if (exposure_id > 0) {
				st3.setLong(1, exposure_id);
				ResultSet rs3 = st3.executeQuery();
				rs3.next();
				mjd_obs = (double) rs3.getDouble(1);
				if (mjd_obs < mjdLo) {mjdLo = mjd_obs;}
				if (mjd_obs > mjdHi) {mjdHi = mjd_obs;}
				rs3.close();
			}

			// Have we encountered this CCD before? If not, add it to the list of parameters 
			//  and find its coordinates on the focal plane (this latter for a QA plot)...
			int found = 0;
			for (iccd = 0; iccd < nccd; iccd++) {
				if (ccdIdArray[iccd] == ccd_number) {
					found = 1;
					break;
				}
			}
			if (found == 0) {
				// if we get here, iccd=nccd
				ccdIdArray[iccd] = ccd_number;
				mStdStarList[iccd] = new ArrayList();
				
				FocalPlaneCoords fpc = this.findFocalPlanePixelCoordinates(ccd_number, image_id, imageType, db);
				nAxis1[iccd] = fpc.getNAxis1();
				nAxis2[iccd] = fpc.getNAxis2();
				xCenter[iccd] = fpc.getXCenter();
				yCenter[iccd] = fpc.getYCenter();
				xAxisMin = Math.min(xAxisMin, xCenter[iccd]-0.5*nAxis1[iccd]);
				xAxisMax = Math.max(xAxisMax, xCenter[iccd]+0.5*nAxis1[iccd]);
				yAxisMin = Math.min(yAxisMin, yCenter[iccd]-0.5*nAxis2[iccd]);
				yAxisMax = Math.max(yAxisMax, yCenter[iccd]+0.5*nAxis2[iccd]);	
				
				nccd++;
			
			}

			// Add mStdStar to the appropriate mStdStarList for the 
			// given CCD
			MatchedStdStarYear1 mStdStar = new MatchedStdStarYear1();
			mStdStar.setStdStarName(stdStarName);
			mStdStar.setFieldName(fieldName);
			mStdStar.setAirmass(airmass);
			mStdStar.setCcd_number(ccd_number);
			mStdStar.setDeltamag(deltamag);
			mStdStar.setDeltamagerr(totMagErr);
			mStdStar.setMjd(mjd_obs);
			mStdStar.setStdmag(stdmag[filterIndex]);
			mStdStar.setStdug(stdmag[0]-stdmag[1]);
			mStdStar.setStdgr(stdmag[1]-stdmag[2]);
			mStdStar.setStdri(stdmag[2]-stdmag[3]);
			mStdStar.setStdiz(stdmag[3]-stdmag[4]);
			mStdStar.setStdzY(stdmag[4]-stdmag[5]);
			mStdStar.setImage_id(image_id);
			mStdStar.setExposure_id(exposure_id);
			mStdStar.setX_image(x_image);
			mStdStar.setY_image(y_image);
			mStdStar.setObject_id(object_id);
			mStdStarList[iccd].add(mStdStar);

			double x = x_image + (xCenter[iccd]-0.5*nAxis1[iccd]);
			double y = y_image + (yCenter[iccd]-0.5*nAxis2[iccd]);

			if (verbose > 1) {
				System.out.println("   " + i + " " + ccd_number + " "
						+ image_id + " " + exposure_id + " " + mjd_obs + " " + 
						+ object_id + " " + x_image + " " + y_image + " " + standard_star_id + " "
						+ stdmag[filterIndex] + " " + stdmag[cFilterIndex] + " " 
						+ instmag + " " + instmagErr + " " + airmass + " " + fieldName); 
			}

			i++;
			
		}

		rs0.close();
		st0.close();
		st1.close();
		st2.close();
		st3.close();

		if (verbose > 2) {
			System.out.println("iccd \t ccdIdArray[iccd]");
			for (iccd = 0; iccd < nccd; iccd++) {
				System.out.println(iccd + "\t" + ccdIdArray[iccd]);
			}
			System.out.println("");
		}

		if (verbose > 0) {
			System.out.println("Fitting data points.");
			System.out.println("");
		}

		if (mjdLo >= 100000) {
			mjdLo = -1;
		}
		if (mjdHi <= 0) {
			mjdHi = -1;
		}
		
		// calculate the number of parameters in the fit...
		// if we are fitting only for the photometric zeropoints
		//  (a_1, ..., a_nccd) and the first-order extinction (k),
		//  the number of free parameters is just nccd plus one.
		int nparam = 2*nccd + 1;
		int nFreeParam = nccd + 1;
		// if we are also fitting for the instrumental color term
		//  ("b term") coefficients (b_1, ..., b_nccd), we have 
		//  an additional nccd free parameters (nFreeParam=nparam).
		if (bsolve) {
			nFreeParam = nFreeParam + nccd;
		}
				
		
		// Set default values for the instrumental color term ("b term") coefficients.
		// If only one value for b has been passed, use this value for all the CCDs;
		//  otherwise, if a list of b values has been passed, use the appropriate value
		//  for each CCD.
		double[] bdefaultValues = new double[nccd];
		double[] bdefaultErrValues = new double[nccd];
		if (colorTermCoeffs.getBccdidArrayList().size() == 1 && 
				Integer.parseInt(colorTermCoeffs.getBccdidArrayList().get(0).toString()) == 0) {
			for (iccd = 0; iccd < nccd; iccd++) {
				bdefaultValues[iccd] = Double.parseDouble(colorTermCoeffs.getBdefaultArrayList().get(0).toString());
				bdefaultErrValues[iccd] = Double.parseDouble(colorTermCoeffs.getBdefaultErrArrayList().get(0).toString());
			}
		} else {
			for (iccd = 0; iccd < nccd; iccd++) {
				int ccdId = ccdIdArray[iccd];
				int index;
				for (index = 0; index < colorTermCoeffs.getBccdidArrayList().size(); index++) {
					int bccdId = Integer.parseInt(colorTermCoeffs.getBccdidArrayList().get(index).toString());
					if (bccdId == ccdId) {
						bdefaultValues[iccd] = Double.parseDouble(colorTermCoeffs.getBdefaultArrayList().get(index).toString());					
						bdefaultErrValues[iccd] = Double.parseDouble(colorTermCoeffs.getBdefaultErrArrayList().get(index).toString());					
					}
				}
			}
		}
		
		
		// initialize a bunch of stuff that sits in the iteration loop
		double[] a = new double[nccd];
		double[] aerr = new double[nccd];
		double[] b = new double[nccd];
		double[] berr = new double[nccd];
		double k = -1;
		double kerr = -1;
		double chi2 = -1;
		double rms = -1;
		int dof = -1;
		int photometricFlag = -1;

		DoubleMatrix2D AA = null;
		DoubleMatrix2D BB = null;
		DoubleMatrix2D XX = null;
		DoubleMatrix2D AAinv = null;
		DoubleMatrix2D II = null;

		// We want to iterate over the solution, culling outliers at each
		// iteration
		for (int iteration = 0; iteration < niterations; iteration++) {

			int iteration1 = iteration + 1;
			if (verbose > 1) {
				System.out.println("   ...  Iteration " + iteration1 + " of "
						+ niterations);
				System.out.println("");
			}

			//(Re-)initialize arrays...
			double[][] array2d = new double[nparam][nparam];
			double[] array1d = new double[nparam];
			
			// Populate arrays containing the matrix to be inverted...
			//   Parameter "0" is k; parameters 1->nccd are a_1, 
			//   ..., a_nccd; parameters nccd+1->2*nccd are b_1, 
			//   ..., b_nccd.
			for (iccd = 0; iccd < nccd; iccd++) {
				int iparam_a = 1 + iccd;
				int iparam_b = 1 + nccd + iccd;
				int size = mStdStarList[iccd].size();
				if (size > 0) {
					for (int j = 0; j < size; j++) {
						
						MatchedStdStarYear1 mStdStar = (MatchedStdStarYear1) mStdStarList[iccd].get(j);
						
						double airmass = mStdStar.getAirmass();
						double deltamag = mStdStar.getDeltamag();
						double error = mStdStar.getDeltamagerr();
						double weight = 1. / (error * error);
						double stdColor = 0.0;
						if (filter.equals("u")) {
							stdColor = mStdStar.getStdug();
						} else if (filter.equals("g") || filter.equals("r")) {
							stdColor = mStdStar.getStdgr();
						} else if (filter.equals("i") || filter.equals("z")) {
							stdColor = mStdStar.getStdiz();
						} else if (filter.equals("Y")) {
							stdColor = mStdStar.getStdzY();
						}

						double deltaStdColor = stdColor-stdColor0;
						
						array2d[0][0] = array2d[0][0] + airmass * airmass * weight;
						array2d[iparam_a][iparam_a] = array2d[iparam_a][iparam_a] + 1 * weight;
						array2d[0][iparam_a] = array2d[0][iparam_a] + airmass * weight;
						array2d[iparam_a][0] = array2d[iparam_a][0] + airmass * weight;
						array1d[0] = array1d[0] + deltamag * airmass * weight;
						array1d[iparam_a] = array1d[iparam_a] + deltamag * weight;

						if (bsolve) {
							array2d[iparam_b][iparam_b] = array2d[iparam_b][iparam_b] 
							                              + deltaStdColor * deltaStdColor * weight;
							array2d[0][iparam_b] = array2d[0][iparam_b] 
							                              + deltaStdColor * airmass * weight;
							array2d[iparam_b][0] = array2d[iparam_b][0] 
							                              + deltaStdColor * airmass * weight;
							array2d[iparam_a][iparam_b] = array2d[iparam_a][iparam_b] 
							                              + deltaStdColor * weight;
							array2d[iparam_b][iparam_a] = array2d[iparam_b][iparam_a] 
							                              + deltaStdColor * weight;
							array1d[iparam_b] = array1d[iparam_b] 
							                              + deltaStdColor * deltamag * weight;
						} else {

							array2d[iparam_b][iparam_b] = 1.0;
							array2d[0][iparam_b] = array2d[0][iparam_b] 
							                               + deltaStdColor * airmass * weight;
							array2d[iparam_b][0] = 0.0;
							array2d[iparam_a][iparam_b] = array2d[iparam_a][iparam_b] 
							                                                + deltaStdColor * weight;
							array2d[iparam_b][iparam_a] = 0.0;
							array1d[iparam_b] = bdefaultValues[iccd];
							
						}
						
					}

				}
				
			}

			if (!ksolve) {
				array2d[0][0] = 1.0;
				array1d[0] = kdefault;
				for (iccd = 0; iccd < nccd; iccd++) {
					int iparam_a = 1 + iccd;
					int iparam_b = 1 + nccd + iccd;
					array2d[0][iparam_a] = 0.0;
					array2d[0][iparam_b] = 0.0;
				}
			}
			
			// Now convert the arrays into matrices.
			AA = new SparseDoubleMatrix2D(nparam, nparam);
			BB = new DenseDoubleMatrix2D(nparam, 1);
			for (int iparam = 0; iparam < nparam; iparam++) {
				BB.set(iparam, 0, array1d[iparam]);
				for (int jparam = 0; jparam < nparam; jparam++) {
					AA.set(iparam, jparam, array2d[iparam][jparam]);
				}
			}

			// Finally, solve the linear equations by matrix inversion...
			if (verbose > 2) {
				System.out.println("Matrix AA: \n" + AA.toString());
				System.out.println("");
			}
			if (verbose > 1) {
				System.out.println("start " + nparam + "x" + nparam
						+ " matrix inversion");
				System.out.println("");
			}
			Algebra alg = new Algebra();
			AAinv = alg.inverse(AA);
			if (verbose > 1) {
				System.out.println(nparam + "x" + nparam
						+ " matrix inversion finished");
				System.out.println("");
			}
			if (verbose > 2) {
				System.out.println("Matrix AAinv: \n" + AAinv.toString());
				System.out.println("");
				II = alg.mult(AA, AAinv);
				System.out.println("Identity Matrix: \n" + II.toString());
				System.out.println("");
			}
			XX = alg.mult(AAinv, BB);
			if (verbose > 2) {
				System.out.println("Matrix XX: \n" + XX.toString());
				System.out.println("");
			}

			// Extract values for k, kerr, a_1, ... , a_nccd, 
			// aerr_1,...,aerr_nccd, b_1, ..., b_nccd, berr_1, 
			// ..., berr_N, from the matrices...
			k = XX.get(0, 0);
			kerr = AAinv.get(0, 0);
			if (ksolve) {
				if (kerr > 0.) {
					kerr = Math.sqrt(kerr);
				} else {
					kerr = -1;
				}
			} else {
				kerr = kdefaultErr;
			}
			for (iccd = 0; iccd < nccd; iccd++) {
				int iparam_a = iccd + 1;
				a[iccd] = XX.get(iparam_a, 0);
				aerr[iccd] = AAinv.get(iparam_a, iparam_a);
				if (aerr[iccd] > 0.) {
					aerr[iccd] = Math.sqrt(aerr[iccd]);
				} else {
					aerr[iccd] = -1;
				}
			}
			for (iccd = 0; iccd < nccd; iccd++) {
				int iparam_b = 1 + nccd + iccd;
				b[iccd] = XX.get(iparam_b, 0);
				if (bsolve) {
					berr[iccd] = AAinv.get(iparam_b, iparam_b);
					if (berr[iccd] > 0.) {
						berr[iccd] = Math.sqrt(berr[iccd]);
					} else {
						berr[iccd] = -1;
					}
				} else  {
					berr[iccd] = bdefaultErrValues[iccd];
				}
			}
			if (verbose > 1) {
				System.out.println("k=  " + k + " +/- " + kerr);
				for (iccd = 0; iccd < nccd; iccd++) {
					int size = mStdStarList[iccd].size();
					System.out.println("a_" + ccdIdArray[iccd] + "= "
							+ a[iccd] + " +/- " + aerr[iccd] + "\t" + "b_"
							+ ccdIdArray[iccd] + "= " + b[iccd] + " +/- "
							+ berr[iccd] + "\t (using " + size
							+ " std stars on this CCD)");
				}
				System.out.println("");
			}

			// Calculate rms and reduced chi2 of solution...
			rms = -1;
			chi2 = -1;
			photometricFlag = -1;
			double sumres2 = 0.;
			double sumchi2 = 0.00;
			int ntot = 0;
			for (iccd = 0; iccd < nccd; iccd++) {
				int size = mStdStarList[iccd].size();
				if (size > 0) {
					ntot = ntot + size;
					for (int j = 0; j < size; j++) {
						MatchedStdStarYear1 mStdStar = (MatchedStdStarYear1) mStdStarList[iccd].get(j);
						double airmass = mStdStar.getAirmass();
						double deltamag = mStdStar.getDeltamag();
						double deltamagerr = mStdStar.getDeltamagerr();
						double stdColor = stdColor0;
						if (filter.equals("u")) {
							stdColor = mStdStar.getStdug();
						} else if (filter.equals("g") || filter.equals("r")) {
							stdColor = mStdStar.getStdgr();
						} else if (filter.equals("i") || filter.equals("z")) {
							stdColor = mStdStar.getStdiz();
						} else if (filter.equals("Y")) {
							stdColor = mStdStar.getStdzY();
						}
						double deltaStdColor = stdColor-stdColor0;
						double res = deltamag - (a[iccd] + k * airmass + b[iccd] * deltaStdColor);
						sumres2 = sumres2 + res * res;
						sumchi2 = sumchi2 + (res / deltamagerr)
								* (res / deltamagerr);
					}
				}
			}
			dof = ntot - nFreeParam;
			if (dof > 0) {
				chi2 = sumchi2 / dof;
				if (sumres2 > 0.) {
					rms = Math.sqrt(sumres2 / dof);
				}
			}
			if (rms > 0. && rms < 0.02) {
				photometricFlag = 1;
			} else {
				photometricFlag = 0;
				System.out.println("QA3BEG rms of solution = " + rms + ", which is considered non-photometric. QA3END");
			}
			if (verbose > 1) {
				System.out.println("        ntot=" + ntot + "  dof=" + dof
						+ "  rms=" + rms + "  chi2=" + chi2);
				System.out.println("");
			}

			// Cull outliers (if this was not the final iteration)
			if (iteration < niterations - 1) {
				// need to work backwards from highest index to lowest...
				if (verbose > 1) {
					System.out.println("        (removing outliers)");
				}
				for (iccd = 0; iccd < nccd; iccd++) {
					int size = mStdStarList[iccd].size();
					if (size > 0) {
						for (int j = 0; j < size; j++) {
							int jj = (size - 1) - j;
							MatchedStdStarYear1 mStdStar = (MatchedStdStarYear1) mStdStarList[iccd].get(jj);
							double airmass = mStdStar.getAirmass();
							double deltamag = mStdStar.getDeltamag();
							double stdColor = 0.0;
							if (filter.equals("u")) {
								stdColor = mStdStar.getStdug();
							} else if (filter.equals("g") || filter.equals("r")) {
								stdColor = mStdStar.getStdgr();
							} else if (filter.equals("i") || filter.equals("z")) {
								stdColor = mStdStar.getStdiz();
							} else if (filter.equals("Y")) {
								stdColor = mStdStar.getStdzY();
							}
							double deltaStdColor = stdColor-stdColor0;
							double res = deltamag - (a[iccd] + k * airmass + b[iccd] * deltaStdColor);
							double resNSigma = res / rms;
							if (Math.abs(res) > nsigma * rms) {
								if (verbose > 1) {
									long object_id = mStdStar.getObject_id();
									System.out
									.println("        Removing outlier (object_id: " 
														 + object_id +  
														 ") on CCD "
														 + ccdIdArray[iccd]
														 + " at airmass "
											             + airmass
											             + " with residual "
											             + res
											             + " (nsigma=" + resNSigma + ")");
								}
								mStdStarList[iccd].remove(jj);
							}
						}
					}
					//if (verbose > 0) {
						//System.out.println("ccd: " + ccdIdArray[iccd] + "\t" + "old size: " + size + "\t new size: " + mStdStarList[iccd].size());
					//}
				}

				System.out.println("");
			}

		}

		if (verbose > 0) {
			System.out.println("");
			System.out.println("Fit completed.");
			System.out.println("");
			System.out.println("Outputting results of fit.");
			System.out.println("");
			System.out.println("Fit Method Name= Matrix Inversion using CERN colt java libraries");
			System.out.println("nite= " + nite);
			System.out.println("MJD range= " + mjdLo + " - " + mjdHi);
			System.out.println("filter= " + filter);
			System.out.println("k=  " + k + " +/- " + kerr);
			for (iccd = 0; iccd < nccd; iccd++) {
				int size = mStdStarList[iccd].size();
				System.out.println("a_" + ccdIdArray[iccd] + "= "
						+ a[iccd] + " +/- " + aerr[iccd] + "\t" + "b_"
						+ ccdIdArray[iccd] + "= " + b[iccd] + " +/- "
						+ berr[iccd] + "\t (using " + size
						+ " std stars on this CCD)");
			}
			System.out.println("rms=" + rms);
			System.out.println("Chi2=" + chi2);
			System.out.println("dof=" + dof);
			System.out.println("");
		}

		
		
		// Prepare output and QA plots... 
		
		if (verbose > 0) {
			System.out.println("Preparing new " + fitTable + " entries and QA plots...");
			System.out.println("");
		}

		if (verbose == 1 && updateDB) {
			System.out.println("Since updateDB option was set, the " + fitTable + " table will be updated.");
			System.out.println("");
		}
		
		// Find latest psmfit_id in database
		Statement st4 = db.createStatement();
		ResultSet rs4 = st4.executeQuery(query4);
		rs4.next();
		double dpsmfit_id = rs4.getDouble(1);
		int psmfit_id = (int) dpsmfit_id;
		if (verbose > 2) {
			System.out.println(dpsmfit_id + " \t " + psmfit_id);
			System.out.println("");
		}
		rs4.close();
		st4.close();

		// there is probably a better way to add a timestamp...
		// java.util.Date d = new java.util.Date();
		Timestamp tt = new Timestamp(date.getTime());

		// We want to use int's rather than booleans for the asolve, 
		// bsolve, and ksolve fields of the PSMFIT table
		int asolveFlag = 0;
		int bsolveFlag = 0;
		int ksolveFlag = 0;
		if (asolve) {asolveFlag = 1;}
		if (bsolve) {bsolveFlag = 1;}
		if (ksolve) {ksolveFlag = 1;}

		
		// These arrays are needed for the creation of the binary FITS table
		//  that will hold the results of the photometric solution.
		//  This is a new feature as of DC5, and is still a little kludgy.
		//  Hopefully, we can clean this up soon.
		int[]    psmfit_idArray       = new int[nccd];
		String[] niteArray            = new String[nccd];
		double[] mjdLoArray           = new double[nccd];
		double[] mjdHiArray           = new double[nccd];
		int[]    ccdId1Array          = new int[nccd];
		double[] kArray               = new double[nccd];
		double[] kerrArray            = new double[nccd];
		double[] rmsArray             = new double[nccd];
		double[] chi2Array            = new double[nccd];
		int[] dofArray                = new int[nccd];
		double[] stdColor0Array       = new double[nccd];
		int[]    photometricFlagArray = new int[nccd];
		String[] psmVersionArray      = new String[nccd];
		String[] timestampArray       = new String[nccd];
		String[] filterArray          = new String[nccd];
		String[] cFilterArray         = new String[nccd];
		String[] runArray             = new String[nccd];
		String[] projectArray         = new String[nccd];
		int[]    asolveFlagArray      = new int[nccd];
		int[]    bsolveFlagArray      = new int[nccd];
		int[]    ksolveFlagArray      = new int[nccd];
		String[] magTypeArray         = new String[nccd];

		
		
		// Loop through the CCDs, preparing the output for each...
		
		for (iccd = 0; iccd < nccd; iccd++) {

			psmfit_id++;

			String values = psmfit_id + ", " + "'" + nite + "', " + mjdLo
					+ ", " + mjdHi + ", " + ccdIdArray[iccd] + ", '"
					+ filter + "', " + a[iccd] + ", " + aerr[iccd] + ", "
					+ b[iccd] + ", " + berr[iccd] + ", " + k + ", " + kerr
					+ ", " + rms + ", " + chi2 + ", " + dof + ", "
					+ photometricFlag + ", '" + psmVersion
					+ "', to_timestamp('" + tt.toString()
					+ "','YYYY-MM-DD HH24:MI:SS.FF3'), '" 
					+ cFilter  + "', " + stdColor0 + ", " 
					+ asolveFlag + ", " + bsolveFlag + ", " + ksolveFlag  + 
					", '" + run + "', '" + project + "', '" + magType + "'";

			// Update the array values for this CCD...
			psmfit_idArray[iccd] = psmfit_id;
			niteArray[iccd] = nite;
			mjdLoArray[iccd] = mjdLo;
			mjdHiArray[iccd] = mjdHi;
			ccdId1Array[iccd] = ccdIdArray[iccd];
			kArray[iccd] = k;
			kerrArray[iccd] = kerr;
			rmsArray[iccd] = rms;
			chi2Array[iccd] = chi2;
			dofArray[iccd] = dof;
			photometricFlagArray[iccd] = photometricFlag;
			psmVersionArray[iccd] = psmVersion;
			filterArray[iccd] = filter;
			cFilterArray[iccd] = cFilter;
			stdColor0Array[iccd] = stdColor0;
			asolveFlagArray[iccd] = asolveFlag;
			bsolveFlagArray[iccd] = bsolveFlag;
			ksolveFlagArray[iccd] = ksolveFlag;
			runArray[iccd] = run;
			projectArray[iccd] = project;
			timestampArray[iccd] = tt.toString();
			magTypeArray[iccd] = magType;


			// Set to true when updating the database table fitTable...
			if (updateDB) {
				if (verbose > 1) {
					System.out.println("Inserting following values into table "
							+ fitTable + " (entry " + psmfit_id + "): ");
					System.out.println(values);
					System.out.println("");
				}
				Statement stmt = db.createStatement();
				stmt.executeUpdate("INSERT INTO " + fitTable + " " + "VALUES ("
						+ values + " )");
				stmt.close();
			} else {
				if (verbose > 1) {
					System.out.println("If the database were being updated, " +
							"the following values would be inserted into table "
							+ fitTable + " (entry " + psmfit_id + "): ");
					System.out.println(values);
					System.out.println("");
				}
			}

			
			// Output a QA plot containing the deltamags vs. airmass for the standard stars
			//  on this CCD...
			if (skipQAPlots == false) {
				
				String qaPlotFile = "PSM_QA_" + nite + filter
				+ ccdIdArray[iccd] + "_" + fitTable + psmfit_id + ".jpg";

				if (verbose > 1) {
					System.out.println("Creating plot " + qaPlotFile + " for CCD "
							+ ccdIdArray[iccd] + "...");
					System.out.println("");
				}

				XYSeries series = new XYSeries("First");

				int size = mStdStarList[iccd].size();
				if (size > 0) {
					for (int j = 0; j < size; j++) {
						MatchedStdStarYear1 mStdStar = (MatchedStdStarYear1) mStdStarList[iccd].get(j);
						double airmass = mStdStar.getAirmass();
						double deltamag = mStdStar.getDeltamag();
						series.add(airmass,deltamag);
					}
				}
				//         Add the series to your data set
				XYSeriesCollection dataset = new XYSeriesCollection();
				dataset.addSeries(series);
				//         Generate the graph
				String title = "Night: " + nite + " Filter: " + filter + " CCD: " + ccdIdArray[iccd];
				JFreeChart chart = ChartFactory.createScatterPlot(
						title, // Title
						"airmass", // x-axis Label
						"instr mag - std mag", // y-axis Label
						dataset, // Dataset
						PlotOrientation.VERTICAL, // Plot Orientation
						false, // Show Legend
						false, // Use tooltips
						false // Configure chart to generate URLs?
				);

				XYPlot plot = (XYPlot) chart.getPlot();
				XYDotRenderer renderer = new XYDotRenderer();
				renderer.setDotHeight(4);
				renderer.setDotWidth(4);
				renderer.setSeriesPaint(0, Color.blue);
				//renderer.setSeriesShape(0, Shape);
				plot.setRenderer(renderer);

				// increase the margins to account for the fact that the auto-range 
				// doesn't take into account the bubble size...
				NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
				domainAxis.setLowerMargin(0.15);
				domainAxis.setUpperMargin(0.15);
				NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
				rangeAxis.setLowerMargin(0.15);
				rangeAxis.setUpperMargin(0.15);

				try {
					ChartUtilities.saveChartAsJPEG(new File(qaPlotFile), chart, 500, 300);
				} catch (IOException e) {
					System.err.println("Problem occurred creating chart.");
				}

			}
	        
		}

		
		// Now actually create the FITS binary table containing the results of the photometric solution...
		

		if (verbose > 0) {
			System.out.println("");
			System.out.println("Creating binary FITS table containing the results of the photometric fit");
			System.out.println("");
		}

		FitsFactory.setUseAsciiTables(false);

		Fits f = new Fits();
		Object[] data  = new Object[]{psmfit_idArray, niteArray, mjdLoArray, mjdHiArray, ccdId1Array,
				filterArray, a, aerr, b, berr, kArray, kerrArray, rmsArray, chi2Array, dofArray, 
				photometricFlagArray, psmVersionArray, timestampArray, cFilterArray, stdColor0Array,
				asolveFlagArray, bsolveFlagArray, ksolveFlagArray, runArray, projectArray, magTypeArray};
		f.addHDU(Fits.makeHDU(data));
		
		BinaryTableHDU bhdu = (BinaryTableHDU) f.getHDU(1);
		bhdu.setColumnName(0,"psmfit_id", null);
		bhdu.setColumnName(1,"nite", null);
		bhdu.setColumnName(2,"mjdLo", "-1 = a dummy value");
		bhdu.setColumnName(3,"mjdHi", "-1 = dummy value");
		bhdu.setColumnName(4,"ccdId", null);
		bhdu.setColumnName(5,"filter", null);
		bhdu.setColumnName(6,"a",null);
		bhdu.setColumnName(7,"aerr",null);
		bhdu.setColumnName(8,"b",null);
		bhdu.setColumnName(9,"berr",null);
		bhdu.setColumnName(10,"k",null);
		bhdu.setColumnName(11,"kerr",null);
		bhdu.setColumnName(12,"rms",null);
		bhdu.setColumnName(13,"chi2","actually, the reduced chi2");
		bhdu.setColumnName(14,"dof",null);
		bhdu.setColumnName(15,"photometricFlag",null);
		bhdu.setColumnName(16,"psmVersion",null);
		bhdu.setColumnName(17,"timestamp",null);
		bhdu.setColumnName(18,"cFilter",null);
		bhdu.setColumnName(19,"stdColor0",null);
		bhdu.setColumnName(20,"asolve",null);
		bhdu.setColumnName(21,"bsolve",null);
		bhdu.setColumnName(22,"ksolve",null);
		bhdu.setColumnName(23,"run",null);
		bhdu.setColumnName(24,"project",null);
		bhdu.setColumnName(25,"magType",null);

		//psmfitInput-$nite-$filter.fits
		String psmfitInputFileName = "psmfitInput_" + nite + filter + ".fits";
		//String psmfitInputFileName = "psmfitInput.fits";
		BufferedFile bf = new BufferedFile(psmfitInputFileName, "rw");
		f.write(bf);
		bf.flush();
		bf.close();

		
		// Close connection to database
		db.close();

		
		// Create and output general QA residual plots and table...
		if (skipQAPlots == false) {

			//Open up the residual table file...
			String resTableFileName = "PSM_QA_res_" + nite + filter + ".txt";
			File resTableFile = new File(resTableFileName);
			FileWriter writer = new FileWriter(resTableFile);
			writer.write("# res \t airmass \t mag \t stdColor \t ccdid \t mjd \t imageid \t exposureid \t X_fp \t Y_fp \t x_ccd \t y_ccd \n");
			
			//Instantiate xy data series for plots...
			XYSeries series1 = new XYSeries("");
			XYSeries series2 = new XYSeries("");
			XYSeries series3 = new XYSeries("");
			XYSeries series4 = new XYSeries("");
			XYSeries series5 = new XYSeries("");
			XYSeries series6 = new XYSeries("");
			//XYSeries series7 = new XYSeries("");
			//XYSeries series8 = new XYSeries("");

			double[] negStarX  = new double[100000];
			double[] negStarY = new double[100000];
			double[] negBubbleSize = new double[100000];
			double[] posStarX  = new double[100000];
			double[] posStarY = new double[100000];
			double[] posBubbleSize = new double[100000];
			int iNegStar = 0;
			int iPosStar = 0;

			for (iccd = 0; iccd < nccd; iccd++) {

				int size = mStdStarList[iccd].size();

				if (size > 0) {

					//series7.add(xCenter[iccd], yCenter[iccd]);

					for (int j = 0; j < size; j++) {

						MatchedStdStarYear1 mStdStar = (MatchedStdStarYear1) mStdStarList[iccd].get(j);

						//String stdStarName = mStdStar.getStdStarName();
						//String fieldName = mStdStar.getFieldName();
						double airmass = mStdStar.getAirmass();
						double mag = mStdStar.getStdmag();
						int ccd_number = mStdStar.getCcd_number();
						double deltamag = mStdStar.getDeltamag();
						double mjd = mStdStar.getMjd();
						long image_id = mStdStar.getImage_id();
						long exposure_id = mStdStar.getExposure_id();

						double x = mStdStar.getX_image() + (xCenter[iccd]-0.5*nAxis1[iccd]);
						double y = mStdStar.getY_image() + (yCenter[iccd]-0.5*nAxis2[iccd]);

						double stdColor = stdColor0;
						if (filter.equals("u")) {
							stdColor = mStdStar.getStdug();
						} else if (filter.equals("g") || filter.equals("r")) {
							stdColor = mStdStar.getStdgr();
						} else if (filter.equals("i") || filter.equals("z")) {
							stdColor = mStdStar.getStdiz();
						} else if (filter.equals("Y")) {
							stdColor = mStdStar.getStdzY();
						}
						double deltaStdColor = stdColor-stdColor0;
						double res = deltamag - (a[iccd] + k * airmass + b[iccd] * deltaStdColor);

						series1.add(airmass,res);
						series2.add(mag,res);
						series3.add(stdColor,res);
						series4.add(ccd_number,res);
						series5.add(mjd,res);
						series6.add(image_id,res);
						if (ccdIdArray[iccd] > -1) {
							//series8.add(x,y);
							// the area of the bubble is proportional to 
							//  the magnitude of the residual...
							double bubbleSize = 0.00;
							if (res != 0.) {
								bubbleSize = 2000.*Math.sqrt(Math.abs(res));
							}
							if (res < 0.) {
								negStarX[iNegStar] = x;
								negStarY[iNegStar] = y;
								negBubbleSize[iNegStar] = bubbleSize;
								iNegStar++;
							} else {
								posStarX[iPosStar] = x;
								posStarY[iPosStar] = y;
								posBubbleSize[iPosStar] = bubbleSize;
								iPosStar++;
							}

						}
						
						//Write entry to residual table file...
						String outputLine = res + "\t" + airmass + "\t" + mag + "\t" + stdColor + "\t" + ccd_number + "\t" + mjd + 
							"\t" + image_id + "\t" + exposure_id + "\t" + x + "\t" + y + "\t" + mStdStar.getX_image() + "\t" + mStdStar.getY_image() + "\n";
						writer.write(outputLine);

					}

				}

			}

			writer.close();

			String qaPlotFile;

			qaPlotFile = "PSM_QA_res_vs_airmass_" + nite + filter + ".jpg";

			if (verbose > 1) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			XYSeriesCollection dataset = new XYSeriesCollection();
			dataset.addSeries(series1);

			String title = "Night: " + nite + " Filter: " + filter;
			JFreeChart chart = ChartFactory.createScatterPlot(
					title, // Title
					"airmass", // x-axis Label
					"mag residual", // y-axis Label
					dataset, // Dataset
					PlotOrientation.VERTICAL, // Plot Orientation
					false, // Show Legend
					false, // Use tooltips
					false // Configure chart to generate URLs?
			);

			XYPlot plot = (XYPlot) chart.getPlot();
			XYDotRenderer renderer = new XYDotRenderer();
			renderer.setDotHeight(4);
			renderer.setDotWidth(4);
			renderer.setSeriesPaint(0, Color.blue);
			//renderer.setSeriesShape(0, Shape);
			plot.setRenderer(renderer);

			// increase the margins to account for the fact that the auto-range 
			// doesn't take into account the bubble size...
			NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
			domainAxis.setLowerMargin(0.15);
			domainAxis.setUpperMargin(0.15);
			NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
			rangeAxis.setLowerMargin(0.15);
			rangeAxis.setUpperMargin(0.15);

			try {
				ChartUtilities.saveChartAsJPEG(new File(qaPlotFile), chart, 500, 300);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}


			qaPlotFile = "PSM_QA_res_vs_mag_" + nite + filter + ".jpg";

			if (verbose > 1) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			dataset = new XYSeriesCollection();
			dataset.addSeries(series2);

			title = "Night: " + nite + " Filter: " + filter;
			chart = ChartFactory.createScatterPlot(
					title, // Title
					"mag", // x-axis Label
					"mag residual", // y-axis Label
					dataset, // Dataset
					PlotOrientation.VERTICAL, // Plot Orientation
					false, // Show Legend
					false, // Use tooltips
					false // Configure chart to generate URLs?
			);

			plot = (XYPlot) chart.getPlot();
			renderer = new XYDotRenderer();
			renderer.setDotHeight(4);
			renderer.setDotWidth(4);
			renderer.setSeriesPaint(0, Color.blue);
			//renderer.setSeriesShape(0, Shape);
			plot.setRenderer(renderer);

			// increase the margins to account for the fact that the auto-range 
			// doesn't take into account the bubble size...
			domainAxis = (NumberAxis) plot.getDomainAxis();
			domainAxis.setLowerMargin(0.15);
			domainAxis.setUpperMargin(0.15);
			rangeAxis = (NumberAxis) plot.getRangeAxis();
			rangeAxis.setLowerMargin(0.15);
			rangeAxis.setUpperMargin(0.15);

			try {
				ChartUtilities.saveChartAsJPEG(new File(qaPlotFile), chart, 500, 300);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}


			qaPlotFile = "PSM_QA_res_vs_color_" + nite + filter + ".jpg";

			if (verbose > 1) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			dataset = new XYSeriesCollection();
			dataset.addSeries(series3);

			title = "Night: " + nite + " Filter: " + filter;
			chart = ChartFactory.createScatterPlot(
					title, // Title
					stdColorName, // x-axis Label
					"mag residual", // y-axis Label
					dataset, // Dataset
					PlotOrientation.VERTICAL, // Plot Orientation
					false, // Show Legend
					false, // Use tooltips
					false // Configure chart to generate URLs?
			);

			plot = (XYPlot) chart.getPlot();
			renderer = new XYDotRenderer();
			renderer.setDotHeight(4);
			renderer.setDotWidth(4);
			renderer.setSeriesPaint(0, Color.blue);
			//renderer.setSeriesShape(0, Shape);
			plot.setRenderer(renderer);

			// increase the margins to account for the fact that the auto-range 
			// doesn't take into account the bubble size...
			domainAxis = (NumberAxis) plot.getDomainAxis();
			domainAxis.setLowerMargin(0.15);
			domainAxis.setUpperMargin(0.15);
			rangeAxis = (NumberAxis) plot.getRangeAxis();
			rangeAxis.setLowerMargin(0.15);
			rangeAxis.setUpperMargin(0.15);

			try {
				ChartUtilities.saveChartAsJPEG(new File(qaPlotFile), chart, 500, 300);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}


			qaPlotFile = "PSM_QA_res_vs_ccd_" + nite + filter + ".jpg";

			if (verbose > 1) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			dataset = new XYSeriesCollection();
			dataset.addSeries(series4);

			title = "Night: " + nite + " Filter: " + filter;
			chart = ChartFactory.createScatterPlot(
					title, // Title
					"CCD number", // x-axis Label
					"mag residual", // y-axis Label
					dataset, // Dataset
					PlotOrientation.VERTICAL, // Plot Orientation
					false, // Show Legend
					false, // Use tooltips
					false // Configure chart to generate URLs?
			);

			plot = (XYPlot) chart.getPlot();
			renderer = new XYDotRenderer();
			renderer.setDotHeight(4);
			renderer.setDotWidth(4);
			renderer.setSeriesPaint(0, Color.blue);
			//renderer.setSeriesShape(0, Shape);
			plot.setRenderer(renderer);

			// increase the margins to account for the fact that the auto-range 
			// doesn't take into account the bubble size...
			domainAxis = (NumberAxis) plot.getDomainAxis();
			domainAxis.setLowerMargin(0.15);
			domainAxis.setUpperMargin(0.15);
			rangeAxis = (NumberAxis) plot.getRangeAxis();
			rangeAxis.setLowerMargin(0.15);
			rangeAxis.setUpperMargin(0.15);

			try {
				ChartUtilities.saveChartAsJPEG(new File(qaPlotFile), chart, 500, 300);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}


			qaPlotFile = "PSM_QA_res_vs_mjd_" + nite + filter + ".jpg";

			if (verbose > 1) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			dataset = new XYSeriesCollection();
			dataset.addSeries(series5);

			title = "Night: " + nite + " Filter: " + filter;
			chart = ChartFactory.createScatterPlot(
					title, // Title
					"mjd_obs", // x-axis Label
					"mag residual", // y-axis Label
					dataset, // Dataset
					PlotOrientation.VERTICAL, // Plot Orientation
					false, // Show Legend
					false, // Use tooltips
					false // Configure chart to generate URLs?
			);

			plot = (XYPlot) chart.getPlot();
			renderer = new XYDotRenderer();
			renderer.setDotHeight(4);
			renderer.setDotWidth(4);
			renderer.setSeriesPaint(0, Color.blue);
			//renderer.setSeriesShape(0, Shape);
			plot.setRenderer(renderer);

			// increase the margins to account for the fact that the auto-range 
			// doesn't take into account the bubble size...
			domainAxis = (NumberAxis) plot.getDomainAxis();
			domainAxis.setLowerMargin(0.15);
			domainAxis.setUpperMargin(0.15);
			rangeAxis = (NumberAxis) plot.getRangeAxis();
			rangeAxis.setLowerMargin(0.15);
			rangeAxis.setUpperMargin(0.15);

			try {
				ChartUtilities.saveChartAsJPEG(new File(qaPlotFile), chart, 500, 300);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}

			qaPlotFile = "PSM_QA_res_vs_imageid_" + nite + filter + ".jpg";

			if (verbose > 1) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			dataset = new XYSeriesCollection();
			dataset.addSeries(series6);

			title = "Night: " + nite + " Filter: " + filter;
			chart = ChartFactory.createScatterPlot(
					title, // Title
					"image_id", // x-axis Label
					"mag residual", // y-axis Label
					dataset, // Dataset
					PlotOrientation.VERTICAL, // Plot Orientation
					false, // Show Legend
					false, // Use tooltips
					false // Configure chart to generate URLs?
			);

			plot = (XYPlot) chart.getPlot();
			renderer = new XYDotRenderer();
			renderer.setDotHeight(4);
			renderer.setDotWidth(4);
			renderer.setSeriesPaint(0, Color.blue);
			//renderer.setSeriesShape(0, Shape);
			plot.setRenderer(renderer);

			// increase the margins to account for the fact that the auto-range 
			// doesn't take into account the bubble size...
			domainAxis = (NumberAxis) plot.getDomainAxis();
			domainAxis.setLowerMargin(0.15);
			domainAxis.setUpperMargin(0.15);
			rangeAxis = (NumberAxis) plot.getRangeAxis();
			rangeAxis.setLowerMargin(0.15);
			rangeAxis.setUpperMargin(0.15);

			try {
				ChartUtilities.saveChartAsJPEG(new File(qaPlotFile), chart, 500, 300);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}


			qaPlotFile = "PSM_QA_res_vs_focalplaneXY_" + nite + filter + ".jpg";

			if (verbose > 1) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			// Define range of the plot
			double axisMin = Math.min(xAxisMin, yAxisMin);
			double axisMax = Math.max(xAxisMax, yAxisMax);
			double extraSpace = 0.10*(axisMax-axisMin);
			axisMin = axisMin - extraSpace;
			axisMax = axisMax + extraSpace;

			// Create dataset.
			DefaultXYZDataset dataxyzset = new DefaultXYZDataset(); 
			double[][] xyzseries0 = new double[][] { negStarX, negStarY, negBubbleSize };
			double[][] xyzseries1 = new double[][] { posStarX, posStarY, posBubbleSize };

			// Create legend to be plotted below the figure.
			double[] xPosLegendList = new double[6]; 
			double[] yPosLegendList = new double[6];
			double[] zPosLegendList = new double[6];
			double[] xNegLegendList = new double[5]; 
			double[] yNegLegendList = new double[5];
			double[] zNegLegendList = new double[5];

			double spacing = (axisMax-axisMin)/12.;
			for (int iii=0; iii<5; iii++) {
				xNegLegendList[iii] = axisMin + (1+iii)*spacing;
				yNegLegendList[iii] = axisMin + 0.5*spacing;
				zNegLegendList[iii] = 2000.*Math.sqrt(Math.abs(2*0.01*(5-iii)));
			}
			for (int iii=0; iii<6; iii++) {
				xPosLegendList[iii] = axisMin + (6+iii)*spacing;
				yPosLegendList[iii] = axisMin + 0.5*spacing;
				if (iii == 0) {
					zPosLegendList[iii] = 0.00;
				} else {
					zPosLegendList[iii] = 2000.*Math.sqrt(Math.abs(2*0.01*(iii)));
				}
			}
			double[][] xyzseries0leg = new double[][] { xNegLegendList, yNegLegendList, zNegLegendList };
			double[][] xyzseries1leg = new double[][] { xPosLegendList, yPosLegendList, zPosLegendList };

			dataxyzset.addSeries("Series 0", xyzseries0);
			dataxyzset.addSeries("Series 1", xyzseries1);
			dataxyzset.addSeries("Series 0 Legend", xyzseries0leg);
			dataxyzset.addSeries("Series 1 Legend", xyzseries1leg);

			//dataset = new XYSeriesCollection();
			//dataset.addSeries(series7);
			//dataset.addSeries(series8);

			title = "Magnitude Residuals vs. Position on Focal Plane\n Night: " + nite + " Filter: " + filter;
			chart = ChartFactory.createBubbleChart(
					title, // Title
					"x [pixels]",   // x-axis Label
					"y [pixels]",  // y-axis Label
					dataxyzset,   // Dataset
					PlotOrientation.VERTICAL, // Plot Orientation
					false,     // Show Legend
					true,      // Use tooltips
					false      // Configure chart to generate URLs?
			);

			plot = (XYPlot) chart.getPlot();

			XYBubbleRenderer xybubblerenderer = new XYBubbleRenderer();

			//modify color of plot symbols 
			XYItemRenderer xyitemRenderer = plot.getRenderer();
			xyitemRenderer.setSeriesPaint(0,Color.red);
			xyitemRenderer.setSeriesPaint(2,Color.red);
			xyitemRenderer.setSeriesPaint(1, Color.green);
			xyitemRenderer.setSeriesPaint(3, Color.green);
			//Paint paint = xyitemRenderer.getSeriesPaint(1);
			//System.out.println("Transparency = " + paint.getTransparency());


			// annotate the plot with the locations and the IDs of the CCDs
			for (iccd=0; iccd<nccd; iccd++) {

				XYBoxAnnotation a1 = new XYBoxAnnotation(
						xCenter[iccd]-0.5*nAxis1[iccd], yCenter[iccd]-0.5*nAxis2[iccd], 
						xCenter[iccd]+0.5*nAxis1[iccd], yCenter[iccd]+0.5*nAxis2[iccd], 
						new BasicStroke(1.5f), Color.black);
				plot.addAnnotation(a1);

				XYTextAnnotation a2 = null;
				//Font font = new Font("SansSerif", Font.PLAIN, 9);
				Font font = new Font("SansSerif", Font.BOLD, 15);
				a2 = new XYTextAnnotation(Integer.toString(ccdIdArray[iccd]), xCenter[iccd], yCenter[iccd]);
				a2.setFont(font);
				a2.setTextAnchor(TextAnchor.CENTER);
				//a2.setTextAnchor(TextAnchor.HALF_ASCENT_LEFT);
				plot.addAnnotation(a2);

			}

			// add annotations for the plot legend
			XYTextAnnotation a3 = null;
			Font font = new Font("SansSerif", Font.PLAIN, 9);

			for (int iii = 0; iii<11; iii++) {
				double resValue = 0.02*(iii-5);
				String resLabel = Double.toString(0.02*(iii-5));
				if (resValue > 0.00) {
					resLabel = "+"+resLabel;
				} 
				//System.out.println(resLabel);
				double xLabel = axisMin + (iii+1)*spacing;
				double yLabel = axisMin + 0.25*spacing;
				a3 = new XYTextAnnotation(resLabel, xLabel, yLabel);
				a3.setFont(font);
				a3.setTextAnchor(TextAnchor.CENTER);
				plot.addAnnotation(a3);

			}

			// increase the margins to account for the fact that the auto-range 
			// doesn't take into account the ccd size...
			domainAxis = (NumberAxis) plot.getDomainAxis();
			domainAxis.setRange(axisMin,axisMax);
			rangeAxis = (NumberAxis) plot.getRangeAxis();
			rangeAxis.setRange(axisMin,axisMax);

			// create and display a frame on screen...
			if (false) {
				ChartFrame frame = new ChartFrame("Focal Plane Plot", chart);
				ChartPanel panel = frame.getChartPanel();
				panel.setFocusable(true);
				panel.setDisplayToolTips(true);
				panel.setPreferredSize(new java.awt.Dimension(680, 680));
				frame.pack();
				frame.setVisible(true);  
				//System.out.println(frame.getChartPanel().getHeight() + " x " + frame.getChartPanel().getWidth());
			}

			//Output plot to JPEG image file...
			try {
				ChartUtilities.saveChartAsJPEG(new File(qaPlotFile), chart, 1200, 1200);
				//ChartUtilities.saveChartAsJPEG(new File(qaPlotFile), chart, 680, 680);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}

		}		
		
		// Finis.
		if (verbose > 0) {
			System.out.println("A binary FITS table containing the results formatted for the " + fitTable + " table can be found in " + psmfitInputFileName);
			System.out.println("");
			System.out.println("That's all, folks!");
			System.out.println("");
		}

	}

	/**
	 * Method findFocalPlanePixelCoordinates tries to find the center and dimensions of 
	 * a CCD on the focal plane measured in pixel coordinates.
	 * (Perhaps this method better belongs in the class FocalPlaneCoords???)
	 * @param int ccd_number, int imageid, String imageType, Connection db
	 * @return Returns the sqlDriver.
	 */
	public FocalPlaneCoords findFocalPlanePixelCoordinates(int ccd_number, long imageid, String imageType, Connection db) throws Exception {

		int nAxis1 = 0;
		int nAxis2 = 0;
		double xCenter = 0.0;
		double yCenter = 0.0;
		
		//With the many new project names, it is easier to assume all but the various BCS projects use the DECam focal plane.
		if (project.equalsIgnoreCase("BCS") || project.equalsIgnoreCase("SCS") || project.equalsIgnoreCase("SPT") || project.equalsIgnoreCase("CPT")) {

			// Assume Blanco+Mosaic2 and use offsets from deprecated wcsoffset table...
			
			if (imageType.equalsIgnoreCase("red")) {
			
				if (ccd_number == 1) {
					
					//CCD 1
					nAxis2  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis1  = (int) Math.round(2*0.0754*3600./0.27);;
					yCenter = -0.1524*3600./0.27;
					xCenter = -0.2344*3600./0.27;

				} else if (ccd_number == 2) {

					//CCD 2
					nAxis2  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis1  = (int) Math.round(2*0.0754*3600./0.27);;
					yCenter = -0.1539*3600./0.27;
					xCenter = -0.0791*3600./0.27;

				} else if (ccd_number == 3) {

					//CCD 3
					nAxis2  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis1  = (int) Math.round(2*0.0754*3600./0.27);;
					yCenter = -0.1545*3600./0.27;
					xCenter = 0.0778*3600./0.27;

				} else if (ccd_number == 4) {

					//CCD 4
					nAxis2  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis1  = (int) Math.round(2*0.0754*3600./0.27);;
					yCenter = -0.1544*3600./0.27;
					xCenter = 0.2327*3600./0.27;

				} else if (ccd_number == 5) {

					//CCD 5
					nAxis2  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis1  = (int) Math.round(2*0.0754*3600./0.27);;
					yCenter = 0.1542*3600./0.27;
					xCenter = -0.2327*3600./0.27;

				} else if (ccd_number == 6) {

					//CCD 6
					nAxis2  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis1  = (int) Math.round(2*0.0754*3600./0.27);;
					yCenter = 0.1545*3600./0.27;
					xCenter = -0.0774*3600./0.27;

				} else if (ccd_number == 7) {

					//CCD 7
					nAxis2  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis1  = (int) Math.round(2*0.0754*3600./0.27);;
					yCenter = 0.1537*3600./0.27;
					xCenter = 0.0774*3600./0.27;

				} else if (ccd_number == 8) {

					//CCD 8
					nAxis2  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis1  = (int) Math.round(2*0.0754*3600./0.27);;
					yCenter = 0.1522*3600./0.27;
					xCenter = 0.2343*3600./0.27;

				}

			} else {

				if (ccd_number == 1) {

					//CCD 1
					nAxis1  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis2  = (int) Math.round(2*0.0754*3600./0.27);;
					xCenter = -0.1524*3600./0.27;
					yCenter = -0.2344*3600./0.27;

				} else if (ccd_number == 2) {

					//CCD 2
					nAxis1  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis2  = (int) Math.round(2*0.0754*3600./0.27);;
					xCenter = -0.1539*3600./0.27;
					yCenter = -0.0791*3600./0.27;

				} else if (ccd_number == 3) {

					//CCD 3
					nAxis1  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis2  = (int) Math.round(2*0.0754*3600./0.27);;
					xCenter = -0.1545*3600./0.27;
					yCenter = 0.0778*3600./0.27;

				} else if (ccd_number == 4) {

					//CCD 4
					nAxis1  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis2  = (int) Math.round(2*0.0754*3600./0.27);;
					xCenter = -0.1544*3600./0.27;
					yCenter = 0.2327*3600./0.27;

				} else if (ccd_number == 5) {

					//CCD 5
					nAxis1  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis2  = (int) Math.round(2*0.0754*3600./0.27);;
					xCenter = 0.1542*3600./0.27;
					yCenter = -0.2327*3600./0.27;

				} else if (ccd_number == 6) {

					//CCD 6
					nAxis1  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis2  = (int) Math.round(2*0.0754*3600./0.27);;
					xCenter = 0.1545*3600./0.27;
					yCenter = -0.0774*3600./0.27;

				} else if (ccd_number == 7) {

					//CCD 7
					nAxis1  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis2  = (int) Math.round(2*0.0754*3600./0.27);;
					xCenter = 0.1537*3600./0.27;
					yCenter = 0.0774*3600./0.27;

				} else if (ccd_number == 8) {

					//CCD 8
					nAxis1  = (int) Math.round(2*0.1508*3600./0.27);
					nAxis2  = (int) Math.round(2*0.0754*3600./0.27);;
					xCenter = 0.1522*3600./0.27;
					yCenter = 0.2343*3600./0.27;

				}

			}
			
		} else {
			
			// Assumes crpix1,crpix2,naxis1,naxis2 are set as in the DECam simulated images
			
			String query = "select crpix1,crpix2,naxis1,naxis2 from image where id = ?";
			PreparedStatement st5 = db.prepareStatement(query);

			st5.setLong(1, imageid);
			ResultSet rs5 = st5.executeQuery();
			rs5.next();
			
			nAxis1  = rs5.getInt("naxis1");
			nAxis2  = rs5.getInt("naxis2");
			xCenter = rs5.getDouble("crpix1") - 0.5*nAxis1;
			yCenter = rs5.getDouble("crpix2") - 0.5*nAxis2;

			rs5.close();
			st5.close();
			
		}
		
		FocalPlaneCoords fpc = new FocalPlaneCoords();
		fpc.setNAxis1(nAxis1);
		fpc.setNAxis2(nAxis2);
		fpc.setXCenter(xCenter);
		fpc.setYCenter(yCenter);
		
		//System.out.println("! ! ! ! ! ! \t" + ccd_number + "\t" + 
		//		nAxis1 + "\t" + nAxis2 + "\t" + 
		//		xCenter + "\t" + yCenter + "\t ! ! ! ! ! !");
		
		return fpc;
		
	}
	
	/**
	 * @return Returns the sqlDriver.
	 */
	public String getSqlDriver() {
		return sqlDriver;
	}

	/**
	 * @param sqlDriver
	 *            The sqlDriver to set.
	 */
	public void setSqlDriver(String sqlDriver) {
		this.sqlDriver = sqlDriver;
	}

	/**
	 * @return Returns the dbName.
	 */
	public String getDbName() {
		return dbName;
	}

	/**
	 * @param dbName
	 *            The dbName to set.
	 */
	public void setDbName(String dbName) {
		this.dbName = dbName;
	}

	/**
	 * @return Returns the project name.
	 */
	public String getProject() {
		return project;
	}

	/**
	 * @param project
	 *            The project name to set.
	 */
	public void setProject(String project) {
		this.project = project;
	}

	/**
	 * @return Returns the url.
	 */
	public String getUrl() {
		return url;
	}

	/**
	 * @param url
	 *            The url to set.
	 */
	public void setUrl(String url) {
		this.url = url;
	}

	/**
	 * @return Returns the user.
	 */
	public String getUser() {
		return user;
	}

	/**
	 * @param user
	 *            The user to set.
	 */
	public void setUser(String user) {
		this.user = user;
	}

	/**
	 * @return Returns the ccdid.
	 */
	public int getCcdid() {
		return ccdid;
	}

	/**
	 * @param ccdid
	 *            The ccdid to set.
	 */
	public void setCcdid(int ccdid) {
		this.ccdid = ccdid;
	}

	/**
	 * @return Returns the mjdHi.
	 */
	public double getMjdHi() {
		return mjdHi;
	}

	/**
	 * @param mjdHi
	 *            The mjdHi to set.
	 */
	public void setMjdHi(double mjdHi) {
		this.mjdHi = mjdHi;
	}

	/**
	 * @return Returns the mjdLo.
	 */
	public double getMjdLo() {
		return mjdLo;
	}

	/**
	 * @param mjdLo
	 *            The mjdLo to set.
	 */
	public void setMjdLo(double mjdLo) {
		this.mjdLo = mjdLo;
	}

	/**
	 * @return Returns the obsTable.
	 */
	public String getObsTable() {
		return obsTable;
	}

	/**
	 * @param obsTable
	 *            The obsTable to set.
	 */
	public void setObsTable(String obsTable) {
		this.obsTable = obsTable;
	}

	/**
	 * @return Returns the imageTable.
	 */
	public String getImageTable() {
		return imageTable;
	}

	/**
	 * @param imageTable
	 *            The imageTable to set.
	 */
	public void setImageTable(String imageTable) {
		this.imageTable = imageTable;
	}

	/**
	 * @return Returns the stdTable.
	 */
	public String getStdTable() {
		return stdTable;
	}

	/**
	 * @param stdTable
	 *            The stdTable to set.
	 */
	public void setStdTable(String stdTable) {
		this.stdTable = stdTable;
	}

	/**
	 * @return Returns the filter.
	 */
	public String getFilter() {
		return filter;
	}

	/**
	 * @param filter
	 *            The filter to set.
	 */
	public void setFilter(String filter) {
		this.filter = filter;
	}

	/**
	 * @return Returns the filterList.
	 */
	public String[] getFilterList() {
		return filterList;
	}

	/**
	 * @param filterList
	 *            The filterList to set.
	 */
	public void setFilterList(String[] filterList) {
		this.filterList = filterList;
	}

	/**
	 * @return Returns the verbose.
	 */
	public int getVerbose() {
		return verbose;
	}

	/**
	 * @param verbose
	 *            The verbose to set.
	 */
	public void setVerbose(int verbose) {
		this.verbose = verbose;
	}

	/**
	 * @return Returns the fitTable.
	 */
	public String getFitTable() {
		return fitTable;
	}

	/**
	 * @param fitTable
	 *            The fitTable to set.
	 */
	public void setFitTable(String fitTable) {
		this.fitTable = fitTable;
	}

	/**
	 * @return Returns the date.
	 */
	public Date getDate() {
		return date;
	}

	/**
	 * @param date
	 *            The date to set.
	 */
	public void setDate(Date date) {
		this.date = date;
	}

	/**
	 * @return Returns the passwd.
	 */
	public String getPasswd() {
		return passwd;
	}

	/**
	 * @param passwd
	 *            The passwd to set.
	 */
	public void setPasswd(String passwd) {
		this.passwd = passwd;
	}

	/**
	 * @return Returns the psmVersion.
	 */
	public String getPsmVersion() {
		return psmVersion;
	}

	/**
	 * @param psmVersion
	 *            The psmVersion to set.
	 */
	public void setPsmVersion(String psmVersion) {
		this.psmVersion = psmVersion;
	}

	/**
	 * @return Returns the magHi.
	 */
	public double getMagHi() {
		return magHi;
	}

	/**
	 * @param magHi
	 *            The magHi to set.
	 */
	public void setMagHi(double magHi) {
		this.magHi = magHi;
	}

	/**
	 * @return Returns the magLo.
	 */
	public double getMagLo() {
		return magLo;
	}

	/**
	 * @param magLo
	 *            The magLo to set.
	 */
	public void setMagLo(double magLo) {
		this.magLo = magLo;
	}

	/**
	 * @return Returns the nite.
	 */
	public String getNite() {
		return nite;
	}

	/**
	 * @param nite
	 *            The nite to set.
	 */
	public void setNite(String nite) {
		this.nite = nite;
	}

	/**
	 * @return Returns the niterations.
	 */
	public int getNiterations() {
		return niterations;
	}

	/**
	 * @param niterations
	 *            The niterations to set.
	 */
	public void setNiterations(int niterations) {
		this.niterations = niterations;
	}

	/**
	 * @return Returns the nsigma.
	 */
	public double getNsigma() {
		return nsigma;
	}

	/**
	 * @param nsigma
	 *            The nsigma to set.
	 */
	public void setNsigma(double nsigma) {
		this.nsigma = nsigma;
	}

	/**
	 * @return Returns the imageNameFilter.
	 */
	public String getImageNameFilter() {
		return imageNameFilter;
	}

	/**
	 * @param imageNameFilter
	 *            The imageNameFilter to set.
	 */
	public void setImageNameFilter(String imageNameFilter) {
		this.imageNameFilter = imageNameFilter;
	}

	/**
	 * @return Returns the imageType.
	 */
	public String getImageType() {
		return imageType;
	}

	/**
	 * @param imageType
	 *            The imageType to set.
	 */
	public void setImageType(String imageType) {
		this.imageType = imageType;
	}

	/**
	 * @return Returns the base, or minimum, error associated with each observed
	 *         mag.
	 */
	public double getBaseMagErr() {
		return baseMagErr;
	}

	/**
	 * @param baseMagErr
	 *            The base, or minimum, error associated with each observed mag,
	 *            to set.
	 */
	public void setBaseMagErr(double baseMagErr) {
		this.baseMagErr = baseMagErr;
	}

	/**
	 * @return Returns the run.
	 */
	public String getRun() {
		return run;
	}

	/**
	 * @param run
	 *            The run to set.
	 */
	public void setRun(String run) {
		this.run = run;
	}

	/**
	 * @return Returns the bsolve.
	 */
	public boolean getBsolve() {
		return bsolve;
	}

	/**
	 * @param bsolve
	 *            The bsolve to set.
	 */
	public void setBsolve(boolean bsolve) {
		this.bsolve = bsolve;
	}

	/**
	 * @return Returns the updateDB.
	 */
	public boolean getUpdateDB() {
		return updateDB;
	}

	/**
	 * @param updateDB
	 *            The updateDB to set.
	 */
	public void setUpdateDB(boolean updateDB) {
		this.updateDB = updateDB;
	}

	/**
	 * @return Returns the ksolve.
	 */
	public boolean getKsolve() {
		return ksolve;
	}

	/**
	 * @param ksolve
	 *            The ksolve to set.
	 */
	public void setKsolve(boolean ksolve) {
		this.ksolve = ksolve;
	}

	/**
	 * @return Returns the stdColor0.
	 */
	public double getStdColor0() {
		return stdColor0;
	}

	/**
	 * @param stdColor0
	 *            The stdColor0 to set.
	 */
	public void setStdColor0(double stdColor0) {
		this.stdColor0 = stdColor0;
	}

	public double getBdefault() {
		return bdefault;
	}

	public void setBdefault(double bdefault) {
		this.bdefault = bdefault;
	}

	public double getBdefaultErr() {
		return bdefaultErr;
	}

	public void setBdefaultErr(double bdefaultErr) {
		this.bdefaultErr = bdefaultErr;
	}

	public double getKdefault() {
		return kdefault;
	}

	public void setKdefault(double kdefault) {
		this.kdefault = kdefault;
	}

	public double getKdefaultErr() {
		return kdefaultErr;
	}

	public void setKdefaultErr(double kdefaultErr) {
		this.kdefaultErr = kdefaultErr;
	}

	/**
	 * @return Returns the colorTermCoeffs.
	 */
	public ColorTermCoeffs getColorTermCoeffs() {
		return colorTermCoeffs;
	}

	/**
	 * @param colorTermCoeffs The colorTermCoeffs to set.
	 */
	public void setColorTermCoeffs(ColorTermCoeffs colorTermCoeffs) {
		this.colorTermCoeffs = colorTermCoeffs;
	}


	public String getMagType() {
		return magType;
	}


	public void setMagType(String magType) {
		this.magType = magType;
	}


	public double getAdefault() {
		return adefault;
	}


	public void setAdefault(double adefault) {
		this.adefault = adefault;
	}


	public double getAdefaultErr() {
		return adefaultErr;
	}


	public void setAdefaultErr(double adefaultErr) {
		this.adefaultErr = adefaultErr;
	}


	public String getImageidExcludeList() {
		return imageidExcludeList;
	}


	public void setImageidExcludeList(String imageidExcludeList) {
		this.imageidExcludeList = imageidExcludeList;
	}

	public boolean getUseOnlyCurrentObjects() {
		return useOnlyCurrentObjects;
	}
	
	public void setUseOnlyCurrentObjects(boolean useOnlyCurrentObjects) {
		this.useOnlyCurrentObjects = useOnlyCurrentObjects;
	}

	public boolean getSkipQAPlots() {
		return skipQAPlots;
	}

	public void setSkipQAPlots(boolean skipQAPlots) {
		this.skipQAPlots = skipQAPlots;
	}

	public int getStandard_set_in() {
		return standard_set_in;
	}

	public void setStandard_set_in(int standard_set_in) {
		this.standard_set_in = standard_set_in;
	}

	public double getStdColorLo() {
		return stdColorLo;
	}

	public void setStdColorLo(double stdColorLo) {
		this.stdColorLo = stdColorLo;
	}

	public double getStdColorHi() {
		return stdColorHi;
	}

	public void setStdColorHi(double stdColorHi) {
		this.stdColorHi = stdColorHi;
	}

}
