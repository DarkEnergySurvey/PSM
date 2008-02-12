package gov.fnal.eag.dtucker.desPhotoStds;

import java.io.IOException;
import java.sql.*;
import java.text.DecimalFormat;

import hep.aida.*;

import java.util.ArrayList;
import java.util.Date;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

// import gov.fnal.eag.dtucker.util.MJD;

/**
 * This class defines methods for solving the photometric zeropoint, a, and the
 * first-order extinction, k, for a given filter for a given night by fitting
 * the following equation for standard star observations: m_inst-m_std = a + kX,
 * where m_inst is the instrumental (observed) mag of a standard star, m_std is
 * the known calibrated mag of the standard star, and X is the airmass of the
 * observation.
 * 
 * @author dtucker
 * 
 */
public class PhotomEqSolver2 {

	// Instance variables dealing with the SQL database
	private String sqlDriver = "oracle.jdbc.driver.OracleDriver";

	private String url = "jdbc:oracle:thin:@//sky.ncsa.uiuc.edu:1521/";

	private String dbName = "des";

	private String user = "dummy";

	private String passwd = "dummy";

	// Instance variables dealing with this run of the Photometric Standards
	// Module
	private Date date = new Date();

	private String psmVersion = null;

	// Instance variables dealing with the observed data to be calibrated
	private String instrument = "DECam";

	private String telescope = "Blanco 4m";

	private String obsTable = "OBJECTS";

	private String filesTable = "FILES";

	private double mjdLo = -1; // no longer selecting on mjd,

	private double mjdHi = -1; // so setting mjdLo=mjdHi=-1

	private int ccdid = 29;

	private String filter = "r";

	private String nite = "2005oct09";

	private double magLo = 17.0;

	private double magHi = 17.1;

	private int niterations = 2; // for outlier removal in fit

	private double nsigma = 5; // for outlier removal in fit

	private String imageType = "REDUCED";

	private String imageNameFilter = "%-0";

	private double baseMagErr = 0.02; // minimum error assoc'd with the
										// observed mags

	private String runiddesc = "%";

	// Instance variables dealing with the standard star list
	private String stdTable = "standard_stars";

	// Instance variables dealing with the output of the fit
	private String fitTable = "psmfit";

	// General instance variables...
	private String[] filterList = { "u", "g", "r", "i", "z" };

	private int verbose = 0;

	public void solve() throws Exception {

		System.out.println("\n\nPhotomEqSolver2");

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
			System.out
					.println("Incompatible filter index.  Throwing Exception!");
			System.out.println("");
			throw new Exception();
		}

		// Stuff from JAS3/JAIDA for making QA plots
		IAnalysisFactory af = IAnalysisFactory.create();
		ITree tree = af.createTreeFactory().create();
		IDataPointSetFactory dpsf = af.createDataPointSetFactory(tree);
		IDataPointSet[] dataPointSet = new IDataPointSet[100];

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
							+ filter + "]  instmag airmass");
		}

		int i = 0;
		double[] stdmag = new double[5];
		double[] stdmagerr = new double[5];

		// Set up an array of parameter IDs
		int paramIdArray[] = new int[100];
		// We already have one parameter, the extinction coefficient k.
		// We will set its parameter ID to -1, to avoid confusion with the CCD
		// zeropoint terms, a1, ..., aN, whose parameter IDs are their
		// CCD_NUMBERs.
		paramIdArray[0] = -1;
		int nparam = 1;
		int iparam = 0;

		Statement st0 = db.createStatement();
		ResultSet rs0 = st0.executeQuery(query0);
		ArrayList mStdStarList = new ArrayList();
		int ccd_numberMin =  10000;
		int ccd_numberMax = -10000;

		while (rs0.next()) {

			int object_id = rs0.getInt("object_id");
			int standard_star_id = rs0.getInt("standard_star_id");

			st1.setInt(1, standard_star_id);
			ResultSet rs1 = st1.executeQuery();
			rs1.next();

			String stdStarName = (String) rs1.getString("name");
			String fieldName = (String) rs1.getString("fieldname");
			//if (fieldName.equalsIgnoreCase("E4_a")) {
				//continue;
			//}
			for (int j = 0; j < filterList.length; j++) {
				String magColumnName = "stdmag_" + filterList[j];
				String magerrColumnName = "stdmagerr_" + filterList[j];
				stdmag[j] = (double) rs1.getFloat(magColumnName);
				stdmagerr[j] = (double) rs1.getFloat(magerrColumnName);
			}

			rs1.close();

			// If the stdmag or the stdmagerr for the filter being solved for
			// doesn't make sense, then skip this std star...
			if (stdmag[filterIndex] < 0. || stdmagerr[filterIndex] < 0) {
				continue;
			}

			double exptime = (double) rs0.getFloat("exptime");
			double instmag = (double) rs0.getFloat("mag_aper_3");
			double zeropoint = (double) rs0.getFloat("zeropoint");
			double airmass = (double) rs0.getFloat("airmass");

			// if (airmass < 1.2) {continue;}

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
			

			// Have we encountered this CCD before? If not, add it to the list
			// of parameters.
			int found = 0;
			for (iparam = 1; iparam < nparam; iparam++) {
				if (paramIdArray[iparam] == ccd_number) {
					found = 1;
					break;
				}
			}
			if (found == 0) {
				paramIdArray[iparam] = ccd_number;
				dataPointSet[iparam] = dpsf.create("dataPointSet" + iparam,
						"CCD " + ccd_number + "\tFilter " + filter + "\tNight "
								+ nite, 2);
				nparam++;
			}

			MatchedStdStar mStdStar = new MatchedStdStar();
			mStdStar.setStdStarName(stdStarName);
			mStdStar.setFieldName(fieldName);
			mStdStar.setAirmass(airmass);
			mStdStar.setCcd_number(ccd_number);
			mStdStar.setDeltamag(deltamag);
			mStdStar.setDeltamagerr(baseMagErr);
			mStdStar.setMjd(0.);
			mStdStar.setStdmag(stdmag[filterIndex]);
			mStdStar.setStdug(stdmag[0]-stdmag[1]);
			mStdStar.setStdgr(stdmag[1]-stdmag[2]);
			mStdStar.setStdri(stdmag[2]-stdmag[3]);
			mStdStar.setStdiz(stdmag[3]-stdmag[4]);
			mStdStarList.add(mStdStar);

			// Add datapoint to the appropriate dataPointSet for the given CCD
			IDataPoint dp = dataPointSet[iparam].addPoint();
			dp.coordinate(0).setValue(airmass);
			dp.coordinate(1).setValue(deltamag);
			dp.coordinate(1).setErrorPlus(baseMagErr);
			dp.coordinate(1).setErrorMinus(baseMagErr);

			if (verbose > 0) {
				System.out.println("   " + i + " " + ccd_number + " "
						+ object_id + " " + standard_star_id + " "
						+ stdmag[filterIndex] + " " + instmag + " " + airmass);
			}

			i++;

		}

		rs0.close();
		st0.close();
		st1.close();

		if (verbose > 1) {
			System.out.println("iparam \t paramIdArray[iparam]");
			for (iparam = 0; iparam < nparam; iparam++) {
				System.out.println(iparam + "\t" + paramIdArray[iparam]);
			}
			System.out.println("");
		}

		if (verbose > 0) {
			System.out.println("Fitting data points.");
			System.out.println("");
		}

		// initialize a bunch of stuff that sits in the iteration loop
		double[] a = new double[100];
		double[] aerr = new double[100];
		double[] b = new double[100];
		double[] berr = new double[100];
		double k = -1;
		double kerr = -1;
		double chi2 = -1;
		double rms = -1;
		int dof = -1;
		int photometricFlag = -1;

		double[][] array2d = new double[100][100];
		double[] array1d = new double[100];

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

			// Populate arrays containing the matrix to be inverted
			for (iparam = 1; iparam < nparam; iparam++) {
				int size = dataPointSet[iparam].size();
				if (size > 0) {
					for (int j = 0; j < size; j++) {
						IDataPoint dp = dataPointSet[iparam].point(j);
						double airmass = dp.coordinate(0).value();
						double deltamag = dp.coordinate(1).value();
						double error = dp.coordinate(1).errorPlus();
						double weight = 1. / (error * error);
						// weight = 1.;
						array2d[0][0] = array2d[0][0] + airmass * airmass
								* weight;
						array2d[iparam][iparam] = array2d[iparam][iparam] + 1
								* weight;
						array2d[0][iparam] = array2d[0][iparam] + airmass
								* weight;
						array2d[iparam][0] = array2d[iparam][0] + airmass
								* weight;
						array1d[0] = array1d[0] + deltamag * airmass * weight;
						array1d[iparam] = array1d[iparam] + deltamag * weight;
					}
				}
			}

			// Now convert the arrays into matrices.
			AA = new SparseDoubleMatrix2D(nparam, nparam);
			BB = new DenseDoubleMatrix2D(nparam, 1);
			for (iparam = 0; iparam < nparam; iparam++) {
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

			// Extract values for k, kerr, a1,...,aN, aerr1,...,aerrN from
			// matrices...
			k = XX.get(0, 0);
			kerr = AAinv.get(0, 0);
			if (kerr > 0.) {
				kerr = Math.sqrt(kerr);
			} else {
				kerr = -1;
			}
			for (iparam = 1; iparam < nparam; iparam++) {
				a[iparam] = XX.get(iparam, 0);
				aerr[iparam] = AAinv.get(iparam, iparam);
				if (aerr[iparam] > 0.) {
					aerr[iparam] = Math.sqrt(aerr[iparam]);
				} else {
					aerr[iparam] = -1;
				}
			}
			if (verbose > 0) {
				System.out.println("k=  " + k + " +/- " + kerr);
				for (iparam = 1; iparam < nparam; iparam++) {
					int size = dataPointSet[iparam].size();
					System.out.println("a_" + paramIdArray[iparam] + "= "
							+ a[iparam] + " +/- " + aerr[iparam] + "\t" + "b_"
							+ paramIdArray[iparam] + "= " + b[iparam] + " +/- "
							+ berr[iparam] + "\t (using " + size
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
			for (iparam = 1; iparam < nparam; iparam++) {
				int size = dataPointSet[iparam].size();
				if (size > 0) {
					ntot = ntot + size;
					for (int j = 0; j < size; j++) {
						IDataPoint dp = dataPointSet[iparam].point(j);
						double airmass = dp.coordinate(0).value();
						double deltamag = dp.coordinate(1).value();
						double res = deltamag - (a[iparam] + k * airmass);
						sumres2 = sumres2 + res * res;
						sumchi2 = sumchi2 + (res / baseMagErr)
								* (res / baseMagErr);
					}
				}
			}
			dof = ntot - nparam;
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
				for (iparam = 1; iparam < nparam; iparam++) {
					int size = dataPointSet[iparam].size();
					if (size > 0) {
						for (int j = 0; j < size; j++) {
							int jj = (size - 1) - j;
							IDataPoint dp = dataPointSet[iparam].point(jj);
							double airmass = dp.coordinate(0).value();
							double deltamag = dp.coordinate(1).value();
							double res = deltamag - (a[iparam] + k * airmass);
							double resNSigma = res / rms;
							if (Math.abs(res) > nsigma * rms) {
								System.out
										.println("        Removing outlier on CCD "
												+ paramIdArray[iparam]
												+ " at airmass "
												+ airmass
												+ " with residual "
												+ res
												+ " (nsigma=" + resNSigma + ")");
								dataPointSet[iparam].removePoint(jj);
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
			System.out
					.println("Fit Method Name= Matrix Inversion using CERN colt java libraries");
			System.out.println("nite= " + nite);
			System.out.println("MJD range= " + mjdLo + " - " + mjdHi);
			System.out.println("filter= " + filter);
			System.out.println("k=  " + k + " +/- " + kerr);
			for (iparam = 1; iparam < nparam; iparam++) {
				int size = dataPointSet[iparam].size();
				System.out.println("a_" + paramIdArray[iparam] + "= "
						+ a[iparam] + " +/- " + aerr[iparam] + "\t" + "b_"
						+ paramIdArray[iparam] + "= " + b[iparam] + " +/- "
						+ berr[iparam] + "\t (using " + size
						+ " std stars on this CCD)");
			}
			System.out.println("rms=" + rms);
			System.out.println("Chi2=" + chi2);
			System.out.println("dof=" + dof);
			System.out.println("");
		}

		if (verbose > 0) {
			System.out
					.println("Preparing to insert new PSMFIT entries into the database and to output QA plots...");
			System.out.println("");
		}

		double value = 0.0;

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

		for (iparam = 1; iparam < nparam; iparam++) {

			psmfit_id++;

			String values = psmfit_id + ", " + "'" + nite + "', " + mjdLo
					+ ", " + mjdHi + ", " + paramIdArray[iparam] + ", '"
					+ filter + "', " + a[iparam] + ", " + aerr[iparam] + ", "
					+ b[iparam] + ", " + berr[iparam] + ", " + k + ", " + kerr
					+ ", " + rms + ", " + chi2 + ", " + dof + ", "
					+ photometricFlag + ", '" + psmVersion
					+ "', to_timestamp('" + tt.toString()
					+ "','YYYY-MM-DD HH24:MI:SS.FF3')";

			if (verbose > 0) {
				System.out.println("Inserting following values into table "
						+ fitTable + " (entry " + psmfit_id + "): ");
				System.out.println(values);
				System.out.println("");
			}

			// Set to false when debugging...
			if (false) {
				Statement stmt = db.createStatement();
				stmt.executeUpdate("INSERT INTO " + fitTable + " " + "VALUES ("
						+ values + " )");
				stmt.close();
			}

			String qaPlotFile = "PSM_QA_" + nite + filter
					+ paramIdArray[iparam] + "_" + fitTable + psmfit_id + ".ps";

			if (verbose > 0) {
				System.out.println("Creating plot " + qaPlotFile + " for CCD "
						+ paramIdArray[iparam] + "...");
				System.out.println("");
			}

			IPlotterFactory plotterFactory = af.createPlotterFactory();
			IPlotter plotter = plotterFactory.create("a plot");
			IPlotterStyle style1 = plotter.region(0).style();
			style1.dataStyle().markerStyle().setParameter("size", "5");
			style1.dataStyle().markerStyle().setParameter("shape", "dot");
			style1.dataStyle().markerStyle().setParameter("color", "blue");
			style1.xAxisStyle().setLabel("airmass");
			style1.yAxisStyle().setLabel("instr mag - std mag");
			style1.setParameter("showStatisticsBox", "false");
			style1.setParameter("showLegend", "false");

			plotter.region(0).setXLimits(0.95, 2.6);
			plotter.region(0).plot(dataPointSet[iparam]);

			IDataPointSet lineDataPointSet = dpsf.create("dataPointSet", "", 2);
			lineDataPointSet.addPoint();
			IDataPoint lineDp = lineDataPointSet.point(0);
			double airmass = 1.0;
			lineDp.coordinate(0).setValue(airmass);
			value = a[iparam] + k * airmass;
			lineDp.coordinate(1).setValue(value);
			lineDataPointSet.addPoint();
			lineDp = lineDataPointSet.point(1);
			airmass = 2.6;
			lineDp.coordinate(0).setValue(airmass);
			value = a[iparam] + k * airmass;
			lineDp.coordinate(1).setValue(value);

			IPlotterStyle style2 = plotterFactory.createPlotterStyle();
			style2.dataStyle().markerStyle().setParameter("size", "1");
			style2.dataStyle().markerStyle().setParameter("shape",
					"horizontalLine");
			style2.dataStyle().markerStyle().setParameter("color", "red");
			// style.dataStyle().setParameter("showDataPoints","false");
			style2.dataStyle().setParameter("connectDataPoints", "true");
			style2.setParameter("showStatisticsBox", "false");
			style2.setParameter("showLegend", "false");

			plotter.region(0).plot(lineDataPointSet, style2);

			// plotter.show();
			try {
				plotter.writeToFile(qaPlotFile, "ps");
			} catch (IOException e) {
				e.printStackTrace();
			}

		}

		// Close connection to database
		db.close();

		if (false) {

			IDataPointSet residDataPointSet = dpsf.create("residDataPointSet",
					nite + " " + filter, 2);
			for (iparam = 1; iparam < nparam; iparam++) {
				int size = dataPointSet[iparam].size();
				if (size > 0) {
					for (int j = 0; j < size; j++) {
						IDataPoint dp = dataPointSet[iparam].point(j);
						double airmass = dp.coordinate(0).value();
						double deltamag = dp.coordinate(1).value();
						double res = deltamag - (a[iparam] + k * airmass);
						IDataPoint resdp = residDataPointSet.addPoint();
						resdp.coordinate(0).setValue(airmass);
						resdp.coordinate(1).setValue(res);
					}
				}
			}

			String qaPlotFile = "PSM_QA_" + nite + filter + ".ps";

			if (verbose > 0) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			IPlotterFactory plotterFactory = af.createPlotterFactory();
			IPlotter plotter = plotterFactory.create("a plot");
			IPlotterStyle style1 = plotter.region(0).style();
			style1.dataStyle().markerStyle().setParameter("size", "5");
			style1.dataStyle().markerStyle().setParameter("shape", "dot");
			style1.dataStyle().markerStyle().setParameter("color", "blue");
			style1.xAxisStyle().setLabel("airmass");
			style1.yAxisStyle().setLabel("mag residuals");
			style1.setParameter("showStatisticsBox", "false");
			style1.setParameter("showLegend", "false");

			plotter.region(0).setXLimits(0.95, 2.6);
			plotter.region(0).plot(residDataPointSet);

			IDataPointSet lineDataPointSet = dpsf.create("dataPointSet", "", 2);
			lineDataPointSet.addPoint();
			IDataPoint lineDp = lineDataPointSet.point(0);
			double airmass = 1.0;
			lineDp.coordinate(0).setValue(airmass);
			value = 0.0;
			lineDp.coordinate(1).setValue(value);
			lineDataPointSet.addPoint();
			lineDp = lineDataPointSet.point(1);
			airmass = 2.6;
			lineDp.coordinate(0).setValue(airmass);
			value = 0.0;
			lineDp.coordinate(1).setValue(value);

			IPlotterStyle style2 = plotterFactory.createPlotterStyle();
			style2.dataStyle().markerStyle().setParameter("size", "1");
			style2.dataStyle().markerStyle().setParameter("shape",
					"horizontalLine");
			style2.dataStyle().markerStyle().setParameter("color", "red");
			// style.dataStyle().setParameter("showDataPoints","false");
			style2.dataStyle().setParameter("connectDataPoints", "true");
			style2.setParameter("showStatisticsBox", "false");
			style2.setParameter("showLegend", "false");

			plotter.region(0).plot(lineDataPointSet, style2);
			// plotter.show();
			try {
				plotter.writeToFile(qaPlotFile, "ps");
			} catch (IOException e) {
				e.printStackTrace();
			}

		} else {

			IDataPointSet residDataPointSet1 = dpsf.create("residuals vs. airmass",
					nite + " " + filter, 2);
			IDataPointSet residDataPointSet2 = dpsf.create("residuals vs. std mag",
					nite + " " + filter, 2);
			IDataPointSet residDataPointSet3 = dpsf.create("residuals vs. std color",
					nite + " " + filter, 2);
			IDataPointSet residDataPointSet4 = dpsf.create("residuals vs. ccd_number",
					nite + " " + filter, 2);
			int size = mStdStarList.size();
			if (size > 0) {
				for (int j = 0; j < size; j++) {
					MatchedStdStar mStdStar = (MatchedStdStar) mStdStarList
							.get(j);
					String stdStarName = mStdStar.getStdStarName();
					String fieldName = mStdStar.getFieldName();
					double airmass = mStdStar.getAirmass();
					double mag = mStdStar.getStdmag();
					double color = mStdStar.getStdgr();
					double ccd_number = mStdStar.getCcd_number();
				    double deltamag = mStdStar.getDeltamag();
					for (iparam = 1; iparam < nparam; iparam++) {
						if (paramIdArray[iparam] == ccd_number) {
							break;
						}
					}

					double res = deltamag - (a[iparam] + k * airmass);
					
					if (airmass < 1.1 && res > 0.5 ) {
						System.out.println(fieldName + "\t" + stdStarName + "\t\t" + ccd_number + "\t" + airmass + "\t" + mag + "\t" + deltamag + "\t" + res);
					}
					
					IDataPoint resdp1 = residDataPointSet1.addPoint();
					resdp1.coordinate(0).setValue(airmass);
					resdp1.coordinate(1).setValue(res);
					
					IDataPoint resdp2 = residDataPointSet2.addPoint();
					resdp2.coordinate(0).setValue(mag);
					resdp2.coordinate(1).setValue(res);
					
					IDataPoint resdp3 = residDataPointSet3.addPoint();
					resdp3.coordinate(0).setValue(color);
					resdp3.coordinate(1).setValue(res);
					
					IDataPoint resdp4 = residDataPointSet4.addPoint();
					resdp4.coordinate(0).setValue(ccd_number);
					resdp4.coordinate(1).setValue(res);
				
				}
			}

			String qaPlotFile;
			
			qaPlotFile = "PSM_QA_res_vs_airmass_" + nite + filter + ".ps";

			if (verbose > 0) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			IPlotterFactory plotterFactory = af.createPlotterFactory();
			IPlotter plotter = plotterFactory.create("a plot");
			IPlotterStyle style1 = plotter.region(0).style();
			style1.dataStyle().markerStyle().setParameter("size", "5");
			style1.dataStyle().markerStyle().setParameter("shape", "dot");
			style1.dataStyle().markerStyle().setParameter("color", "blue");
			style1.xAxisStyle().setLabel("airmass");
			style1.yAxisStyle().setLabel("mag residuals");
			style1.setParameter("showStatisticsBox", "false");
			style1.setParameter("showLegend", "false");

			plotter.region(0).setXLimits(0.95, 2.6);
			//plotter.region(0).setYLimits(-0.5, 0.5);
			plotter.region(0).plot(residDataPointSet1);

			IDataPointSet lineDataPointSet = dpsf.create("dataPointSet", "", 2);
			lineDataPointSet.addPoint();
			IDataPoint lineDp = lineDataPointSet.point(0);
			lineDp.coordinate(0).setValue(1.0);
			lineDp.coordinate(1).setValue(0.0);
			lineDataPointSet.addPoint();
			lineDp = lineDataPointSet.point(1);
			lineDp.coordinate(0).setValue(2.6);
			lineDp.coordinate(1).setValue(0.0);

			IPlotterStyle style2 = plotterFactory.createPlotterStyle();
			style2.dataStyle().markerStyle().setParameter("size", "1");
			style2.dataStyle().markerStyle().setParameter("shape",
					"horizontalLine");
			style2.dataStyle().markerStyle().setParameter("color", "red");
			// style.dataStyle().setParameter("showDataPoints","false");
			style2.dataStyle().setParameter("connectDataPoints", "true");
			style2.setParameter("showStatisticsBox", "false");
			style2.setParameter("showLegend", "false");

			plotter.region(0).plot(lineDataPointSet, style2);
			// plotter.show();
			try {
				plotter.writeToFile(qaPlotFile, "ps");
			} catch (IOException e) {
				e.printStackTrace();
			}


			qaPlotFile = "PSM_QA_res_vs_mag_" + nite + filter + ".ps";

			if (verbose > 0) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			plotterFactory = af.createPlotterFactory();
			plotter = plotterFactory.create("a plot");
			style1 = plotter.region(0).style();
			style1.dataStyle().markerStyle().setParameter("size", "5");
			style1.dataStyle().markerStyle().setParameter("shape", "dot");
			style1.dataStyle().markerStyle().setParameter("color", "blue");
			style1.xAxisStyle().setLabel("std mag ("+filter+")");
			style1.yAxisStyle().setLabel("mag residuals");
			style1.setParameter("showStatisticsBox", "false");
			style1.setParameter("showLegend", "false");

			plotter.region(0).setXLimits(14.0, 20.0);
			//plotter.region(0).setYLimits(-0.5, 0.5);
			plotter.region(0).plot(residDataPointSet2);

			lineDataPointSet = dpsf.create("dataPointSet", "", 2);
			lineDataPointSet.addPoint();
			lineDp = lineDataPointSet.point(0);
			lineDp.coordinate(0).setValue(14.0);
			lineDp.coordinate(1).setValue(0.0);
			lineDataPointSet.addPoint();
			lineDp = lineDataPointSet.point(1);
			lineDp.coordinate(0).setValue(20.0);
			lineDp.coordinate(1).setValue(0.0);

			style2 = plotterFactory.createPlotterStyle();
			style2.dataStyle().markerStyle().setParameter("size", "1");
			style2.dataStyle().markerStyle().setParameter("shape",
					"horizontalLine");
			style2.dataStyle().markerStyle().setParameter("color", "red");
			// style.dataStyle().setParameter("showDataPoints","false");
			style2.dataStyle().setParameter("connectDataPoints", "true");
			style2.setParameter("showStatisticsBox", "false");
			style2.setParameter("showLegend", "false");

			plotter.region(0).plot(lineDataPointSet, style2);
			// plotter.show();
			try {
				plotter.writeToFile(qaPlotFile, "ps");
			} catch (IOException e) {
				e.printStackTrace();
			}


			qaPlotFile = "PSM_QA_res_vs_gr_" + nite + filter + ".ps";

			if (verbose > 0) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			plotterFactory = af.createPlotterFactory();
			plotter = plotterFactory.create("a plot");
			style1 = plotter.region(0).style();
			style1.dataStyle().markerStyle().setParameter("size", "5");
			style1.dataStyle().markerStyle().setParameter("shape", "dot");
			style1.dataStyle().markerStyle().setParameter("color", "blue");
			style1.xAxisStyle().setLabel("std color (g-r)");
			style1.yAxisStyle().setLabel("mag residuals");
			style1.setParameter("showStatisticsBox", "false");
			style1.setParameter("showLegend", "false");

			plotter.region(0).setXLimits(-0.5, 1.8);
			//plotter.region(0).setYLimits(-0.5, 0.5);
			plotter.region(0).plot(residDataPointSet3);

			lineDataPointSet = dpsf.create("dataPointSet", "", 2);
			lineDataPointSet.addPoint();
			lineDp = lineDataPointSet.point(0);
			lineDp.coordinate(0).setValue(-0.5);
			lineDp.coordinate(1).setValue(0.0);
			lineDataPointSet.addPoint();
			lineDp = lineDataPointSet.point(1);
			lineDp.coordinate(0).setValue(1.8);
			lineDp.coordinate(1).setValue(0.0);

			style2 = plotterFactory.createPlotterStyle();
			style2.dataStyle().markerStyle().setParameter("size", "1");
			style2.dataStyle().markerStyle().setParameter("shape",
					"horizontalLine");
			style2.dataStyle().markerStyle().setParameter("color", "red");
			// style.dataStyle().setParameter("showDataPoints","false");
			style2.dataStyle().setParameter("connectDataPoints", "true");
			style2.setParameter("showStatisticsBox", "false");
			style2.setParameter("showLegend", "false");

			plotter.region(0).plot(lineDataPointSet, style2);
			// plotter.show();
			try {
				plotter.writeToFile(qaPlotFile, "ps");
			} catch (IOException e) {
				e.printStackTrace();
			}


			qaPlotFile = "PSM_QA_res_vs_ccd_" + nite + filter + ".ps";

			if (verbose > 0) {
				System.out.println("Creating plot " + qaPlotFile + "...");
				System.out.println("");
			}

			plotterFactory = af.createPlotterFactory();
			plotter = plotterFactory.create("a plot");
			style1 = plotter.region(0).style();
			style1.dataStyle().markerStyle().setParameter("size", "5");
			style1.dataStyle().markerStyle().setParameter("shape", "dot");
			style1.dataStyle().markerStyle().setParameter("color", "blue");
			style1.xAxisStyle().setLabel("ccd_number");
			style1.yAxisStyle().setLabel("mag residuals");
			style1.setParameter("showStatisticsBox", "false");
			style1.setParameter("showLegend", "false");

			plotter.region(0).setXLimits(ccd_numberMin, ccd_numberMax);
			//plotter.region(0).setYLimits(-0.5, 0.5);
			plotter.region(0).plot(residDataPointSet4);

			lineDataPointSet = dpsf.create("dataPointSet", "", 2);
			lineDataPointSet.addPoint();
			lineDp = lineDataPointSet.point(0);
			lineDp.coordinate(0).setValue(0.);
			lineDp.coordinate(1).setValue(0.0);
			lineDataPointSet.addPoint();
			lineDp = lineDataPointSet.point(1);
			lineDp.coordinate(0).setValue(65.);
			lineDp.coordinate(1).setValue(0.0);

			style2 = plotterFactory.createPlotterStyle();
			style2.dataStyle().markerStyle().setParameter("size", "1");
			style2.dataStyle().markerStyle().setParameter("shape",
					"horizontalLine");
			style2.dataStyle().markerStyle().setParameter("color", "red");
			// style.dataStyle().setParameter("showDataPoints","false");
			style2.dataStyle().setParameter("connectDataPoints", "true");
			style2.setParameter("showStatisticsBox", "false");
			style2.setParameter("showLegend", "false");

			plotter.region(0).plot(lineDataPointSet, style2);
			// plotter.show();
			try {
				plotter.writeToFile(qaPlotFile, "ps");
			} catch (IOException e) {
				e.printStackTrace();
			}


		}

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

}
