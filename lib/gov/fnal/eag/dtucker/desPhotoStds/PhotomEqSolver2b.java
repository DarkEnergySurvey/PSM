package gov.fnal.eag.dtucker.desPhotoStds;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.sql.*;
import java.util.ArrayList;
import java.util.Date;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYDotRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;


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
public class PhotomEqSolver2b {

	// Instance variables dealing with the SQL database
	private String sqlDriver = "oracle.jdbc.driver.OracleDriver";
	private String url = "jdbc:oracle:thin:@charon.ncsa.uiuc.edu:1521:";
	private String dbName = "des";
	private String user = "dummy";
	private String passwd = "dummy";
	private String obsTable = "OBJECTS";
	private String filesTable = "FILES";
	private String stdTable = "standard_stars";
	private String fitTable = "psmfit";

	// Instance variables dealing with this run of the Photometric 
	// Standards Module
	private String[] filterList = { "u", "g", "r", "i", "z" };
	private int verbose = 0;
	private String psmVersion = "v2b";
	private Date date = new Date();

	// Instance variables dealing with the observed data to be 
	// calibrated
	private String instrument = "DECam";
	private String telescope = "Blanco 4m";
	private double mjdLo = -1; // no longer selecting on mjd,
	private double mjdHi = -1; // so setting mjdLo=mjdHi=-1
	private int ccdid = 0;     // ccdid=0 means all ccds
	private String filter = "r";
	private String nite = "2005oct09";
	private double magLo = 15.0;
	private double magHi = 18.0;
	private String imageType = "remap";
	private String imageNameFilter = "%";
	private String runiddesc = "%";

	// Instance variables dealing with the fit
	private int niterations = 3; // for outlier removal in fit
	private double nsigma = 2.5; // for outlier removal in fit
	private boolean debug = false; //set to false if updating fitTable
	
	private double stdColor0   = 0.00;
	private double baseMagErr = 0.02; // minimum error assoc'd with the
                                      // observed mag
	
	private boolean asolve      = true; // we don't currently use
	private double  adefault    = 0.00; // these a term private 
	private double  adefaultErr = 0.01; // variables except for PSMFIT

	private boolean bsolve      = false;
	private double  bdefault    = 0.00;
	private double  bdefaultErr = 0.01;

	private boolean ksolve      = false;
	private double  kdefault    = 0.18;
	private double  kdefaultErr = 0.02;
	
	
	public void solve() throws Exception {

		System.out.println("\n\nPhotomEqSolver2b");

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
		if (verbose > 0) {
			System.out.println("filterIndex for filter " + filter + " is "
					+ filterIndex + ".");
			System.out.println("");
		}
		if (filterIndex < 0) {
			System.out.println("Incompatible filter index.  Throwing Exception!");
			System.out.println("");
			throw new Exception();
		}

		//Parameters needed for determination of b terms...
		String cFilter = "";
		String stdColorName = "";
		if (filter.equals("g")) {
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
		}	
		System.out.println("filter is " + filter + " cFilter is " + cFilter);
		int cFilterIndex = -1;
		for (int j = 0; j < filterList.length; j++) {
			if (cFilter.equals(filterList[j])) {
				cFilterIndex = j;
				break;
			}
		}
		if (verbose > 0) {
			System.out.println("cFilterIndex for cFilter " + cFilter + " is "
					+ cFilterIndex + ".");
			System.out.println("");
		}
		if (cFilterIndex < 0) {
			System.out.println("Incompatible cFilter index.  Throwing Exception!");
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

		String query0 = "SELECT * FROM table(fPhotoMatch(" + "'" + imageType
				+ "', '" + imageNameFilter + "', " + "'" + nite + "', '"
				+ filter + "', " + ccdid + ", " + magLo + ", " + magHi + ", '"
				+ instrument + "', '" + telescope + "', '" + runiddesc + "'))";
		if (verbose > 0) {
			System.out.println("query0 = " + query0);
			System.out.println("");
		}

		String query1 = "SELECT * FROM " + stdTable
				+ " WHERE standard_star_id = ?";
		PreparedStatement st1 = db.prepareStatement(query1);
		if (verbose > 0) {
			System.out.println("query1 = " + query1);
			System.out.println("");
		}

		String query2 = "SELECT CCD_NUMBER FROM " + filesTable + " f, "
				+ obsTable + " o WHERE o.imageid=f.imageid AND o.object_id = ?";
		PreparedStatement st2 = db.prepareStatement(query2);
		if (verbose > 0) {
			System.out.println("query2 = " + query2);
			System.out.println("");
		}

		if (verbose > 0) {
			System.out
					.println("Point  ccd_number object_id  standard_star_id  stdmag["
							+ filter + "]  instmag airmass  fieldName");
		}

		int i = 0;
		double[] stdmag = new double[5];
		double[] stdmagerr = new double[5];

		// Set up an array of ccd IDs...
		int ccdIdArray[] = new int[100];
		// iccd is the index of the ccdIdArray, nccd is the total 
		// number of ccds in the fit...
		int iccd = 0;
		int nccd = 0;
		// ccd_numberMin and ccd_numberMax are the minimum and maximum 
		// values of CCD_NUMBER for those ccds in the fit (these values 
		// are used by one of the QA plots)
		int ccd_numberMin =  10000;
		int ccd_numberMax = -10000;

		Statement st0 = db.createStatement();
		ResultSet rs0 = st0.executeQuery(query0);
		ArrayList[] mStdStarList = new ArrayList[100];

		while (rs0.next()) {

			int object_id = rs0.getInt("object_id");
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
			if (stdmag[filterIndex] < 0. || stdmagerr[filterIndex] < 0.) {
				continue;
			}
			// If the stdmag or the stdmagerr for the cFilter for 
			// the filter being solved for doesn't make sense, then 
			// skip this std star...
			if (stdmag[cFilterIndex] < 0. || stdmagerr[cFilterIndex] < 0.) {
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
			double instmag = (double) rs0.getFloat("mag_aper_3");
			double zeropoint = (double) rs0.getFloat("zeropoint");
			double airmass = (double) rs0.getFloat("airmass");

			//if (airmass < 1.05) {continue;}
			
			instmag = instmag - zeropoint;
			if (exptime > 0.) {
				instmag = instmag + 2.5 * 0.4342944819 * Math.log(exptime);
			}
			double deltamag = instmag - stdmag[filterIndex];

			// Find on which CCD this star lies...
			st2.setInt(1, object_id);
			ResultSet rs2 = st2.executeQuery();
			rs2.next();
			int ccd_number = (int) rs2.getInt("ccd_number");
			if (ccd_number < ccd_numberMin) {ccd_numberMin = ccd_number;}
			if (ccd_number > ccd_numberMax) {ccd_numberMax = ccd_number;}
			

			// Have we encountered this CCD before? If not, add it to 
			// the list of parameters.
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
				nccd++;
			}

			// Add mStdStar to the appropriate mStdStarList for the 
			// given CCD
			MatchedStdStar mStdStar = new MatchedStdStar();
			mStdStar.setStdStarName(stdStarName);
			mStdStar.setFieldName(fieldName);
			mStdStar.setAirmass(airmass);
			mStdStar.setCcd_number(ccd_number);
			mStdStar.setDeltamag(deltamag);
			mStdStar.setDeltamagerr(totMagErr);
			mStdStar.setMjd(0.);
			mStdStar.setStdmag(stdmag[filterIndex]);
			mStdStar.setStdug(stdmag[0]-stdmag[1]);
			mStdStar.setStdgr(stdmag[1]-stdmag[2]);
			mStdStar.setStdri(stdmag[2]-stdmag[3]);
			mStdStar.setStdiz(stdmag[3]-stdmag[4]);
			mStdStarList[iccd].add(mStdStar);

			if (verbose > 0) {
				System.out.println("   " + i + " " + ccd_number + " "
						+ object_id + " " + standard_star_id + " "
						+ stdmag[filterIndex] + " " + instmag + " " + airmass + " " + fieldName);
			}

			i++;
			
		}

		rs0.close();
		st0.close();
		st1.close();

		if (verbose > 1) {
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

		// calculate the number of parameters in the fit...
		// if we are fitting only for the photometric zeropoints
		//  (a_1, ..., a_nccd) and the first-order extinction (k),
		//  the number of parameters is just nccd plus one.
		int nparam = 2*nccd + 1;
		int nFreeParam = nccd + 1;
		if (bsolve) {
			nFreeParam = nFreeParam + nccd;
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

		//double[][] array2d = new double[nparam][nparam];
		//double[] array1d = new double[nparam];

		DoubleMatrix2D AA = null;
		DoubleMatrix2D BB = null;
		DoubleMatrix2D XX = null;
		DoubleMatrix2D AAinv = null;
		DoubleMatrix2D II = null;

		// We want to iterate over the solution, culling outliers at each
		// iteration
		for (int iteration = 0; iteration < niterations; iteration++) {

			int iteration1 = iteration + 1;
			if (verbose > 0) {
				System.out.println("   ...  Iteration " + iteration1 + " of "
						+ niterations);
				System.out.println("");
			}

			//(Re-)nitialize arrays...
			double[][] array2d = new double[nparam][nparam];
			double[] array1d = new double[nparam];
			
			// Populate arrays containing the matrix to be inverted...
			//   Parameter "0" is k; parameters 1->nparam are a_1, 
			//   ..., a_nccd.
			for (iccd = 0; iccd < nccd; iccd++) {
				int iparam_a = 1 + iccd;
				int iparam_b = 1 + nccd + iccd;
				int size = mStdStarList[iccd].size();
				if (size > 0) {
					for (int j = 0; j < size; j++) {
						
						MatchedStdStar mStdStar = (MatchedStdStar) mStdStarList[iccd].get(j);
						
						double airmass = mStdStar.getAirmass();
						double deltamag = mStdStar.getDeltamag();
						double error = mStdStar.getDeltamagerr();
						double weight = 1. / (error * error);
						double stdColor = 0.0;
						if (filter.equals("g") || filter.equals("r")) {
							stdColor = mStdStar.getStdgr();
						} else if (filter.equals("i") || filter.equals("z")) {
							stdColor = mStdStar.getStdiz();
						}
						double deltaStdColor = stdColor-stdColor0;
						
						// weight = 1.;
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
							array1d[iparam_b] = bdefault;
							
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

			// Finally, solve the linear equations...
			if (verbose > 1) {
				System.out.println("Matrix AA: \n" + AA.toString());
				System.out.println("");
			}
			if (verbose > 0) {
				System.out.println("start " + nparam + "x" + nparam
						+ " matrix inversion");
				System.out.println("");
			}
			Algebra alg = new Algebra();
			AAinv = alg.inverse(AA);
			if (verbose > 0) {
				System.out.println(nparam + "x" + nparam
						+ " matrix inversion finished");
				System.out.println("");
			}
			if (verbose > 1) {
				System.out.println("Matrix AAinv: \n" + AAinv.toString());
				System.out.println("");
				II = alg.mult(AA, AAinv);
				System.out.println("Identity Matrix: \n" + II.toString());
				System.out.println("");
			}
			XX = alg.mult(AAinv, BB);
			if (verbose > 1) {
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
					berr[iccd] = bdefaultErr;
				}
			}
			if (verbose > 0) {
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
						MatchedStdStar mStdStar = (MatchedStdStar) mStdStarList[iccd].get(j);
						double airmass = mStdStar.getAirmass();
						double deltamag = mStdStar.getDeltamag();
						double deltamagerr = mStdStar.getDeltamagerr();
						double stdColor = stdColor0;
						if (filter.equals("g") || filter.equals("r")) {
							stdColor = mStdStar.getStdgr();
						} else if (filter.equals("i") || filter.equals("z")) {
							stdColor = mStdStar.getStdiz();
						}
						double deltaStdColor = stdColor-stdColor0;
						double res = deltamag - (a[iccd] + k * airmass + b[iccd] * deltaStdColor);
						sumres2 = sumres2 + res * res;
						sumchi2 = sumchi2 + (res / deltamagerr)
								* (res / deltamagerr);
					}
				}
			}
			//dof = ntot - nparam;
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
			}
			if (verbose > 0) {
				System.out.println("        ntot=" + ntot + "  dof=" + dof
						+ "  rms=" + rms + "  chi2=" + chi2);
				System.out.println("");
			}

			// cull outliers (if this was not the final iteration)
			if (iteration < niterations - 1) {
				// need to work backwards from highest index to lowest...
				if (verbose > 0) {
					System.out.println("        (removing outliers)");
				}
				for (iccd = 0; iccd < nccd; iccd++) {
					int size = mStdStarList[iccd].size();
					if (size > 0) {
						for (int j = 0; j < size; j++) {
							int jj = (size - 1) - j;
							MatchedStdStar mStdStar = (MatchedStdStar) mStdStarList[iccd].get(jj);
							double airmass = mStdStar.getAirmass();
							double deltamag = mStdStar.getDeltamag();
							//double stdColor = mStdStar.getStdgr();
							double stdColor = 0.0;
							if (filter.equals("g") || filter.equals("r")) {
								stdColor = mStdStar.getStdgr();
							} else if (filter.equals("i") || filter.equals("z")) {
								stdColor = mStdStar.getStdiz();
							}
							double deltaStdColor = stdColor-stdColor0;
							double res = deltamag - (a[iccd] + k * airmass + b[iccd] * deltaStdColor);
							double resNSigma = res / rms;
							if (Math.abs(res) > nsigma * rms) {
								System.out
										.println("        Removing outlier on CCD "
												+ ccdIdArray[iccd]
												+ " at airmass "
												+ airmass
												+ " with residual "
												+ res
												+ " (nsigma=" + resNSigma + ")");
								mStdStarList[iccd].remove(jj);
							}
						}
					}
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

		if (verbose > 0) {
			System.out.println("Preparing new " + fitTable + " entries and QA plots...");
			System.out.println("");
		}


		// Find latest psmfit_id in database
		String query3 = "SELECT max(psmfit_id) FROM " + fitTable;
		if (verbose > 0) {
			System.out.println("query3 = " + query3);
			System.out.println("");
		}
		Statement st3 = db.createStatement();
		ResultSet rs3 = st3.executeQuery(query3);
		rs3.next();
		double dpsmfit_id = rs3.getDouble(1);
		int psmfit_id = (int) dpsmfit_id;
		if (verbose > 1) {
			System.out.println(dpsmfit_id + " \t " + psmfit_id);
			System.out.println("");
		}

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
		
		
		for (iccd = 0; iccd < nccd; iccd++) {

			psmfit_id++;

			String values = psmfit_id + ", " + "'" + nite + "', " + mjdLo
					+ ", " + mjdHi + ", " + ccdIdArray[iccd] + ", '"
					+ filter + "', " + a[iccd] + ", " + aerr[iccd] + ", "
					+ b[iccd] + ", " + berr[iccd] + ", " + k + ", " + kerr
					+ ", " + rms + ", " + chi2 + ", " + dof + ", "
					+ photometricFlag + ", '" + psmVersion
					+ "', to_timestamp('" + tt.toString()
					+ "','YYYY-MM-DD HH24:MI:SS.FF3'), '" + cFilter  + "', " 
					+ stdColor0 + ", " + asolveFlag + ", " + bsolveFlag + ", " + ksolveFlag;

			//String values = psmfit_id + ", " + "'" + nite + "', " + mjdLo
			//		+ ", " + mjdHi + ", " + ccdIdArray[iccd] + ", '"
			//		+ filter + "', " + a[iccd] + ", " + aerr[iccd] + ", "
			//		+ b[iccd] + ", " + berr[iccd] + ", " + k + ", " + kerr
			//		+ ", " + rms + ", " + chi2 + ", " + dof + ", "
			//		+ photometricFlag + ", '" + psmVersion
			//		+ "', to_timestamp('" + tt.toString()
			//		+ "','YYYY-MM-DD HH24:MI:SS.FF3')";


			// Set to false when debugging...
			if (!debug) {
				if (verbose > 0) {
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
				if (verbose > 0) {
					System.out.println("If the database were being updated, " +
							"the following values would be inserted into table "
							+ fitTable + " (entry " + psmfit_id + "): ");
					System.out.println(values);
					System.out.println("");
				}
			}

			String qaPlotFile = "PSM_QA_" + nite + filter
					+ ccdIdArray[iccd] + "_" + fitTable + psmfit_id + ".jpg";

			if (verbose > 0) {
				System.out.println("Creating plot " + qaPlotFile + " for CCD "
						+ ccdIdArray[iccd] + "...");
				System.out.println("");
			}

			XYSeries series = new XYSeries("First");
			
			int size = mStdStarList[iccd].size();
			if (size > 0) {
				for (int j = 0; j < size; j++) {
					MatchedStdStar mStdStar = (MatchedStdStar) mStdStarList[iccd].get(j);
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

		// Close connection to database
		db.close();

		
		// Create and output general QA residual plots...
		
		XYSeries series1 = new XYSeries("");
		XYSeries series2 = new XYSeries("");
		XYSeries series3 = new XYSeries("");
		XYSeries series4 = new XYSeries("");
		
		
		for (iccd = 0; iccd < nccd; iccd++) {
			
			int size = mStdStarList[iccd].size();
			
			if (size > 0) {
				
				for (int j = 0; j < size; j++) {
					
					MatchedStdStar mStdStar = (MatchedStdStar) mStdStarList[iccd].get(j);
					
					String stdStarName = mStdStar.getStdStarName();
					String fieldName = mStdStar.getFieldName();
					double airmass = mStdStar.getAirmass();
					double mag = mStdStar.getStdmag();
					double ccd_number = mStdStar.getCcd_number();
					double deltamag = mStdStar.getDeltamag();
					for (iccd = 0; iccd < nccd; iccd++) {
						if (ccdIdArray[iccd] == ccd_number) {
							break;
						}
					}
					double stdColor = stdColor0;
					if (filter.equals("g") || filter.equals("r")) {
						stdColor = mStdStar.getStdgr();
					} else if (filter.equals("i") || filter.equals("z")) {
						stdColor = mStdStar.getStdiz();
					}
					double deltaStdColor = stdColor-stdColor0;
					double res = deltamag - (a[iccd] + k * airmass + b[iccd] * deltaStdColor);
					
					//if (airmass < 1.1 && res > 0.5 ) {
						//System.out.println(fieldName + "\t" + stdStarName + "\t\t" + ccd_number + "\t" + airmass + "\t" + mag + "\t" + deltamag + "\t" + res);
					//}
					
					series1.add(airmass,res);
					series2.add(mag,res);
					series3.add(stdColor,res);
					series4.add(ccd_number,res);
					
				}
				
			}
			
		}
		
		String qaPlotFile;
		
		qaPlotFile = "PSM_QA_res_vs_airmass_" + nite + filter + ".jpg";
		
		if (verbose > 0) {
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
		
		if (verbose > 0) {
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
		
		if (verbose > 0) {
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
		
		if (verbose > 0) {
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
		
		
		// Finis.
		if (verbose > 0) {
			System.out.println("That's all, folks!");
			System.out.println("");
		}

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
	 * @return Returns the instrument.
	 */
	public String getInstrument() {
		return instrument;
	}

	/**
	 * @param instrument
	 *            The instrument to set.
	 */
	public void setInstrument(String instrument) {
		this.instrument = instrument;
	}

	/**
	 * @return Returns the telescope.
	 */
	public String getTelescope() {
		return telescope;
	}

	/**
	 * @param telescope
	 *            The telescope to set.
	 */
	public void setTelescope(String telescope) {
		this.telescope = telescope;
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
	 * @return Returns the filesTable.
	 */
	public String getFilesTable() {
		return filesTable;
	}

	/**
	 * @param filesTable
	 *            The filesTable to set.
	 */
	public void setFilesTable(String filesTable) {
		this.filesTable = filesTable;
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
	 * @return Returns the runiddesc.
	 */
	public String getRuniddesc() {
		return runiddesc;
	}

	/**
	 * @param runiddesc
	 *            The runiddesc to set.
	 */
	public void setRuniddesc(String runiddesc) {
		this.runiddesc = runiddesc;
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
	 * @return Returns the debug.
	 */
	public boolean getDebug() {
		return debug;
	}

	/**
	 * @param debug
	 *            The debug to set.
	 */
	public void setDebug(boolean debug) {
		this.debug = debug;
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

}
