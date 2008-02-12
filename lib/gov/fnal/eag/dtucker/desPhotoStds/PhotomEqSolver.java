package gov.fnal.eag.dtucker.desPhotoStds;

import java.io.IOException;
import java.sql.*;
import java.text.DecimalFormat;

import hep.aida.*;
import java.util.Date;
//import gov.fnal.eag.dtucker.util.MJD;


/**
 * This class defines methods for solving the photometric zeropoint, a, 
 * and the first-order extinction, k, for a given filter for a given night
 * by fitting the following equation for standard star observations:
 *    m_inst-m_std = a + kX,
 * where m_inst is the instrumental (observed) mag of a standard star,
 * m_std is the known calibrated mag of the standard star, and 
 * X is the airmass of the observation.
 * @author dtucker
 *
 */
public class PhotomEqSolver {
    
    //Instance variables dealing with the SQL database
    private String sqlDriver = "oracle.jdbc.driver.OracleDriver";
    private String url       = "jdbc:oracle:thin:@//sky.ncsa.uiuc.edu:1521/";
    private String dbName    = "des";
    private String user      = "dummy";
    private String passwd    = "dummy";
  
    //Instance variables dealing with this run of the Photometric Standards Module
    private Date date = new Date();
    private String psmVersion = null;
    
    //Instance variables dealing with the observed data to be calibrated
    private String instrument = "DECam";
    private String telescope  = "Blanco 4m";
    private String obsTable   = "OBJECTS";
    private String filesTable = "FILES";
    private double mjdLo      = -1;     //no longer selecting on mjd, 
    private double mjdHi      = -1;     //so setting mjdLo=mjdHi=-1
    private int    ccdid      = 29;
    private String filter     = "r";
    private String nite       = "2005oct09";
    private double magLo      = 17.0;
    private double magHi      = 17.1;
    private int niterations   = 2;      //for outlier removal in fit
    private double nsigma     = 5;    //for outlier removal in fit
    private String imageType  = "REDUCED";
    private String imageNameFilter = "%-0";
   
    //Instance variables dealing with the standard star list
    private String stdTable  = "des_stripe82_stds_v1";
    
    //Instance variables dealing with the output of the fit
    private String fitTable  = "psmfit";
    
    //General instance variables...
    private String[] filterList = {"u", "g", "r", "i", "z"};
    private int verbose = 0;
    
    public void solve () throws Exception {
        
        if (verbose > 0) {
            System.out.println("");
            System.out.println("The beginning...");
            System.out.println("");
        }
            
        int filterIndex = -1;
        for (int j=0; j < filterList.length; j++) {
            if (filter.equals(filterList[j])) {
                filterIndex = j;
                break;
            }
        }
        if (verbose > 0) {
            System.out.println("filterIndex for filter " + filter + " is " + filterIndex + ".");
            System.out.println("");
        }
        if (filterIndex < 0) {
            System.out.println("Incompatible filter index.  Throwing Exception!");
            System.out.println("");
            throw new Exception();
        }
        
        
        //Stuff from JAS3/JAIDA for doing fits and for making QA plots
        IAnalysisFactory af = IAnalysisFactory.create();
        ITree tree = af.createTreeFactory().create();
        IDataPointSetFactory dpsf = af.createDataPointSetFactory(tree);
        IFunctionFactory funcF = af.createFunctionFactory(tree);
        IFitFactory fitF = af.createFitFactory();
        IFitter fitter = fitF.createFitter("Chi2", "uncmin");
        IDataPointSet dataPointSet = dpsf.create("dataPointSet", "two dimensional IDataPointSet", 2);
        
        //Establish connection to database
        if (verbose > 0) {
            System.out.println("Establishing connection to database.");
            System.out.println("");
        }
        Class.forName(sqlDriver);
        String fullURL = url+dbName;
        Connection db = DriverManager.getConnection(fullURL, user, passwd);
                
        String query0 = "SELECT * FROM table(fPhotoMatch("+
        		"'"+imageType+"', '"+imageNameFilter+"', " +
        		"'"+nite+"', '"+filter+"', "+ccdid+", "+magLo+", "+magHi + 
        		", '"+instrument+"', '"+telescope+"'))";
        if (verbose > 0) {  
            System.out.println("query0 = " + query0);
            System.out.println("");
        }

        String query1 = "SELECT * FROM " + stdTable + " WHERE stripe82_id = ?";
        PreparedStatement st1 = db.prepareStatement(query1);
        if (verbose > 0) {
        	System.out.println("query1 = " + query1);
            System.out.println("");
        }
        
       if (verbose > 0) {
    	   System.out.println("Point  object_id  stripe82_id  stdmag["+filter+"]  instmag airmass");
       }
       
        Statement st0 = db.createStatement();
        ResultSet rs0 = st0.executeQuery(query0);

        int i = 0;
        double[] stdmag    = new double[5];
        double[] stdmagerr = new double[5];

        while (rs0.next()) {

            int object_id   = rs0.getInt("object_id");
            int stripe82_id = rs0.getInt("stripe82_id");

            st1.setInt(1,stripe82_id);
            ResultSet rs1 = st1.executeQuery();
            rs1.next();

            for (int j=0; j<filterList.length; j++) {
                String magColumnName    = "stdmag_"+filterList[j];
                String magerrColumnName = "stdmagerr_"+filterList[j];
                stdmag[j]    = (double) rs1.getFloat(magColumnName);
                stdmagerr[j] = (double) rs1.getFloat(magerrColumnName);
            }

            rs1.close();
            
            // If the stdmag or the stdmagerr for the filter being solved for
            // doesn't make sense, then skip this std star...
            if (stdmag[filterIndex] < 0. || stdmagerr[filterIndex] < 0) {
                continue;
            }
            
            double exptime   = (double) rs0.getFloat("exptime");
            double instmag   = (double) rs0.getFloat("mag_aper_3");
            double zeropoint = (double) rs0.getFloat("zeropoint");
            double airmass   = (double) rs0.getFloat("airmass");
            
            instmag = instmag - zeropoint;
            if (exptime > 0.) {
               instmag = instmag + 2.5*0.4342944819*Math.log(exptime);
            }
            double deltamag = instmag-stdmag[filterIndex];
            
            dataPointSet.addPoint();
            IDataPoint dp = dataPointSet.point(i);
            dp.coordinate(0).setValue(airmass);
            dp.coordinate(1).setValue(deltamag);
            dp.coordinate(1).setErrorPlus(0.02);
            dp.coordinate(1).setErrorMinus(0.02);
            if (verbose > 0) {
                System.out.println("   " + i + 
                			":  " + object_id + " " + stripe82_id + 
                            " " + stdmag[filterIndex] + 
                            " " + instmag +
                            " " + airmass);
            }

            i++;
            
        }
        
        rs0.close();
        st0.close();
        st1.close();

        if (verbose > 0) {
            System.out.println("");
            System.out.println("Plotting data points");
            System.out.println("");
        }
        
        IPlotter plotter = af.createPlotterFactory().create("a plot");
        plotter.region(0).plot(dataPointSet);
        //plotter.show();
        
        if (verbose > 0) {
            System.out.println("Fitting data points.");
            System.out.println("");
        }
        
        // initialize a bunch of stull that sits in the iteration loop
        IFunction line = null;
        IFitResult result = null;
        double a = -1;
        double aerr = -1; 
        double b = -1;
        double berr = -1;
        double k = -1;
        double kerr = -1;
        double chi2 = -1;
        double rms = -1;
        int dof  = -1;
        int photometricFlag = -1;
        
        for (int iteration=0; iteration < niterations; iteration++) {
            
            int iteration1 = iteration + 1;
            if (verbose > 0) {
                System.out.println("   ...  Iteration " + iteration1 + " of " + niterations);
            }
            
            line = funcF.createFunctionByName("line", "p1");
            result = fitter.fit(dataPointSet, line);
            
            String names[] = result.fittedParameterNames();
            a = result.fittedParameter(names[0]);
            b = 0.0;
            k = result.fittedParameter(names[1]);
            
            double [] errors = result.errors();
            aerr = errors[0];
            berr = 0.0;
            kerr = errors[1];
            
            chi2 = result.quality();
            rms = -1.0;
            dof = result.ndf();
            photometricFlag = -1;
            
            //Calculate rms of solution...
            //(probably should do calculation with dof, but use ntot for now)
            double sum  = 0.;
            double sum2 = 0.;
            int ntot = dataPointSet.size();
            if (ntot > 0) {
                for (int j=0; j < ntot; j++) {
                    IDataPoint dp = dataPointSet.point(j);
                    double airmass  = dp.coordinate(0).value();
                    double deltamag = dp.coordinate(1).value();
                    double res = deltamag - (a + k*airmass);
                    sum  = sum  + res;
                    sum2 = sum2 + res*res;
                }
                double ave = sum/ntot;
                double ave2 = sum2/ntot;
                rms = ave2 - ave*ave;
                if (rms > 0.) {
                	rms = Math.sqrt(rms);
                } else {
                	rms = 99.;
                }
                if (rms < 0.02) {
                    photometricFlag = 1;
                } else {
                    photometricFlag = 0;
                }
                if (verbose > 0) {
                	System.out.println("        ntot=" + ntot + "  rms=" + rms);
                }
                
                //cull outliers (if this was not the final iteration)
                if (iteration < niterations-1) {
                    //need to work backwards from highest index to lowest... 
                    if (verbose > 0) {
                        System.out.println("        (removing outliers)");
                    }
                    for (int j=0; j < ntot; j++) {
                        int jj = (ntot-1)-j;
                        IDataPoint dp = dataPointSet.point(jj);
                        double airmass  = dp.coordinate(0).value();
                        double deltamag = dp.coordinate(1).value();
                        double res = deltamag - (a + k*airmass);
                        double resNSigma = res/rms;
                        if (Math.abs(res) > nsigma*rms) {
                        	System.out.println("        Removing outlier at airmass " + airmass + 
                        			" with residual " + res + " (nsigma=" + resNSigma + ")" );
                            dataPointSet.removePoint(jj);
                        }
                    }                    
                }
                
            }

        }
        
        if (verbose > 0) {
            System.out.println("");
            System.out.println("Fit completed.");
            System.out.println("");
            System.out.println("Outputting results of fit.");
            System.out.println("");
            System.out.println("Fit Method Name= " + result.fitMethodName());
            System.out.println("nite= " + nite);
            System.out.println("MJD range= " + mjdLo + " - " + mjdHi);
            System.out.println("filter= " + filter);
            System.out.println("a= " + a + ", aerr= " + aerr);
            System.out.println("b= " + b + ", berr= " + berr);
            System.out.println("k= " + k + ", kerr= " + kerr);
            System.out.println("rms=" + rms);
            System.out.println("Chi2=" + chi2);
            System.out.println("dof=" + dof);
            System.out.println("");
        }
        
        //Find latest psmfit_id
        String query2 = "SELECT max(psmfit_id) FROM " + fitTable;
        if (verbose > 0) {  
            System.out.println("query2 = " + query0);
            System.out.println("");
        }
        Statement st2 = db.createStatement();
        ResultSet rs2 = st2.executeQuery(query2);
        rs2.next();
        double dpsmfit_id = rs2.getDouble(1);
        int psmfit_id = (int) dpsmfit_id;
        System.out.println(dpsmfit_id + " /t " + psmfit_id);
        //int psmfit_id = rs2.getInt(1);
        psmfit_id = psmfit_id + 1;
        
        // there is probably a better way to add a timestamp...
        //java.util.Date d = new java.util.Date();
        Timestamp tt = new Timestamp(date.getTime());
        
        String values = psmfit_id + ", " + 
        				"'" + nite + "', " + mjdLo  + ", " + mjdHi + ", " + 
                        ccdid + ", '" + filter + "', " +
                        a + ", " + aerr + ", " + 
                        b + ", " + berr + ", " + 
                        k + ", " + kerr + ", " + 
                        rms + ", " + chi2 + ", " + dof + ", " + 
                        photometricFlag + ", '" + psmVersion + 
                        "', to_timestamp('" + tt.toString() + "','YYYY-MM-DD HH24:MI:SS.FF3')";
        
        if (verbose > 0) {
            System.out.println("Inserting following values into table " + fitTable + " (entry " + 
            		psmfit_id + "): ");
            System.out.println(values);
            System.out.println("");
        }
        
        Statement stmt = db.createStatement();
        stmt.executeUpdate("INSERT INTO " + fitTable + " " + 
                "VALUES (" + values + " )");
        stmt.close();
        
        //Close connection to database
        db.close();
        
        
        String qaPlotFile = "PSM_QA_"+nite+filter+ccdid+"_"+fitTable+psmfit_id+".ps";
        if (verbose > 0) {
        	System.out.println("Plotting fit to file " + qaPlotFile + ".");
        	System.out.println("");
        }
        plotter.region(0).plot(result.fittedFunction());
        try {
            plotter.writeToFile(qaPlotFile, "ps");
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
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
     * @param sqlDriver The sqlDriver to set.
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
     * @param dbName The dbName to set.
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
     * @param instrument The instrument to set.
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
     * @param telescope The telescope to set.
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
     * @param url The url to set.
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
     * @param user The user to set.
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
     * @param ccdid The ccdid to set.
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
     * @param mjdHi The mjdHi to set.
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
     * @param mjdLo The mjdLo to set.
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
     * @param obsTable The obsTable to set.
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
     * @param filesTable The filesTable to set.
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
     * @param stdTable The stdTable to set.
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
     * @param filter The filter to set.
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
     * @param filterList The filterList to set.
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
     * @param verbose The verbose to set.
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
     * @param fitTable The fitTable to set.
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
     * @param date The date to set.
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
     * @param passwd The passwd to set.
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
     * @param psmVersion The psmVersion to set.
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
     * @param magHi The magHi to set.
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
     * @param magLo The magLo to set.
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
     * @param nite The nite to set.
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
     * @param niterations The niterations to set.
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
     * @param nsigma The nsigma to set.
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
     * @param imageNameFilter The imageNameFilter to set.
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
     * @param imageType The imageType to set.
     */
	public void setImageType(String imageType) {
		this.imageType = imageType;
	}
	
}
