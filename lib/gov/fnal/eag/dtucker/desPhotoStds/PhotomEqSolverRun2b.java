/**
 * PhotomEqSolverRun2b.java   The DES Collaboration, 17 July  2007
 */

package gov.fnal.eag.dtucker.desPhotoStds;

import gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolver2b;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Date;
import java.util.StringTokenizer;

import jargs.gnu.CmdLineParser;


/**
 * @author dtucker
 *
 */
public class PhotomEqSolverRun2b {
	
//	private String[] parameterList = {
//			"url",
//			"dbName",
//			"user",
//			"passwd",
//			"instrument",
//			"telescope",
//			"nite",
//			"filter",
//			"stdColor0",
//			"ccdid",
//			"magLo",
//			"magHi",
//			"niterations",
//			"nsigma",
//			"imageType",
//			"imageNameFilter",
//			"runiddesc",
//			"psmVersion",
//			"bsolve",
//			"bdefault",
//			"bdefaultErr",
//			"ksolve",
//			"kdefault",
//			"kdefaultErr",
//			"debug",
//			"verbose"
//	};
	
	private static void printUsage(boolean error) {

		PhotomEqSolver2b ph = new PhotomEqSolver2b();
        
		// Grab default values from PhotomEqSolver2b...
		String urlDefault             = ph.getUrl();
        String dbNameDefault          = ph.getDbName();
        String userDefault            = ph.getUser();
        String passwdDefault          = ph.getPasswd();
        String instrumentDefault      = ph.getInstrument();
        String telescopeDefault       = ph.getTelescope();
        String niteDefault            = ph.getNite();
        String filterDefault          = ph.getFilter();
        double stdColor0Default       = ph.getStdColor0();
        int ccdidDefault              = ph.getCcdid();
        double magLoDefault           = ph.getMagLo();
        double magHiDefault           = ph.getMagHi();
        int niterationsDefault        = ph.getNiterations();
        double nsigmaDefault          = ph.getNsigma();
        String imageTypeDefault       = ph.getImageType();
        String imageNameFilterDefault = ph.getImageNameFilter();
        String runiddescDefault       = ph.getRuniddesc();
        String psmVersionDefault      = ph.getPsmVersion();
        boolean bsolveDefault         = ph.getBsolve();
        double bdefaultDefault        = ph.getBdefault();
        double bdefaultErrDefault     = ph.getBdefaultErr();
        boolean ksolveDefault         = ph.getKsolve();
        double kdefaultDefault        = ph.getKdefault();
        double kdefaultErrDefault     = ph.getKdefaultErr();
        boolean debugDefault          = ph.getDebug();
        int verboseDefault            = ph.getVerbose();
		
        // Create output message...
        String message = 
       			"\nUsage:  java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRun2b [OPTIONS] \n\n" +
				"   OPTIONS: \n" + 
				"   --paramFile VALUE             name of (optional) parameter file \n" + 
				"                                 values for parameters in the paramFile will override values hardwired into the PSM2b code, \n" +
				"                                 but values passed as optional parameters in the command line will override the values in \n" +
				"                                 the paramFile \n" +
				"   --url VALUE                   URL of database                       [default: " + urlDefault + "] \n" + 
				"   --dbName VALUE                name of database                      [default: " + dbNameDefault + "] \n" + 
				"   -u VALUE, --user VALUE        username for db access                [default: " + userDefault + "] \n" + 
				"   -p VALUE, --passwd VALUE      password for db access                [default: " + passwdDefault + "] \n" + 
				"   -i VALUE, --instrument VALUE  name of instrument in db              [default: " + instrumentDefault + "] \n" +
				"   -t VALUE, --telescope VALUE   name of telescope in db               [default: " + telescopeDefault + "] \n" +
				"   -n VALUE, --nite VALUE        name of night in db                   [default: " + niteDefault + "] \n" +
				"   -f VALUE, --filter VALUE      name of filter in db                  [default: " + filterDefault + "] \n" +
				"   -c VALUE, --ccdid VALUE       ccd number (0=all ccds)               [default: " + ccdidDefault + "] \n" + 
				"   --stdColor0 VALUE             zeropoint color in photometric eqn.   [default: " + stdColor0Default + "] \n" +  
				"   --magLo VALUE                 mag limit (bright)                    [default: " + magLoDefault + "] \n" +
				"   --magHi VALUE                 mag limit (faint)                     [default: " + magHiDefault + "] \n" +
				"   --niter VALUE                 # of interations to the fit           [default: " + niterationsDefault + "] \n" + 
				"   --nsigma VALUE                # of sigma for outlier rejection      [default: " + nsigmaDefault + "] \n" +
				"   --imageType VALUE             image type for std star fields        [default: " + imageTypeDefault + "] \n" +
				"   --imageNameFilter VALUE       image name filter for std star fields [default: " + imageNameFilterDefault + "] \n" +
				"   --runiddesc VALUE             runiddesc for std star fields         [default: " + runiddescDefault + "] \n" + 
				"   --psmVersion VALUE            version of PSM                        [default: " + psmVersionDefault + "] \n" + 	
				"   --bsolve                      include this flag to solve for instrumental color (b) terms \n" +
				"   --bdefault VALUE              default value for b                   [default: " + bdefaultDefault + "] \n" +
				"   --bdefaultErr VALUE           1sigma error in default value for b   [default: " + bdefaultErrDefault + "] \n" +
				"   --ksolve                      include this flag to solve for the first-order extinction (k) term \n" +	
				"   --kdefault VALUE              default value for k                   [default: " + kdefaultDefault + "] \n" +
				"   --kdefaultErr VALUE           1sigma error in default value for k   [default: " + kdefaultErrDefault + "] \n" +
				"   -v VALUE, --verbose VALUE     verbosity level (0, 1, 2, ...)        [default: " + verboseDefault + "] \n" + 
				"   -d, --debug                   include this flag if the database is not to be updated \n" +
				"   -h, --help                    this message \n\n" + 
				"   Example 1: \n" +
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRun2b --url jdbc:oracle:thin:@charon.ncsa.uiuc.edu:1521: --dbName des  -u myUserName -p myPassword -i Mosaic2 -t \"Blanco 4m\" -n bcs061223 -f g --ccdid 0 --magLo 15.0 --magHi 18.0 --niter 3 --nsigma 2.5 --imageType remap --imageNameFilter % --runiddesc % --psmVersion v2b --debug -v 2 --bsolve --ksolve \n\n" +
				"   Example 2: \n" +
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRun2b --user myUserName -p myPassword -i Mosaic2 -n bcs061223 -f g --debug -v 2 --ksolve \n\n" + 
				"   Example 3: \n" +
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRun2b --user myUserName -p myPassword --paramFile /home/myname/psmParam.par \n\n" + 
				"   Example 4: \n" +
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRun2b --user myUserName -p myPassword -i Mosaic2 -n bcs061223 -f g --debug -v 2 --ksolve  --paramFile /home/myname/psmParam.par \n\n" + 
				"   Example 5: \n" + 
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRun2b --help \n\n";


				
        // Output usage message...
        if (error) {
        	// use System.err.println if this is part of an error message...
        	System.err.println(message);
        } else {
        	// use System.out.println otherwise (e.g., if -h or --help were invoked)...
        	System.out.println(message);
        }

	}

    public static void main (String[] args) throws Exception {

    	System.out.println("PhotomEqSolverRun2b");
    	
        System.out.print("arglist:  ");
        for (int i=0; i< args.length; i++) {
            System.out.print(args[i] + " ");
        }
        System.out.print("\n");
        
        // Instantiate an instance of the PhotomEqSolver2b class...
        PhotomEqSolver2b ph = new PhotomEqSolver2b();

        // Grab default values for parameters from PhotomEqSolver2b...
		String urlDefault             = ph.getUrl();
        String dbNameDefault          = ph.getDbName();
        String userDefault            = ph.getUser();
        String passwdDefault          = ph.getPasswd();
        String instrumentDefault      = ph.getInstrument();
        String telescopeDefault       = ph.getTelescope();
        String niteDefault            = ph.getNite();
        String filterDefault          = ph.getFilter();
        double stdColor0Default       = ph.getStdColor0();
        int ccdidDefault              = ph.getCcdid();
        double magLoDefault           = ph.getMagLo();
        double magHiDefault           = ph.getMagHi();
        int niterationsDefault        = ph.getNiterations();
        double nsigmaDefault          = ph.getNsigma();
        String imageTypeDefault       = ph.getImageType();
        String imageNameFilterDefault = ph.getImageNameFilter();
        String runiddescDefault       = ph.getRuniddesc();
        String psmVersionDefault      = ph.getPsmVersion();
        boolean bsolveDefault         = ph.getBsolve();
        double bdefaultDefault        = ph.getBdefault();
        double bdefaultErrDefault     = ph.getBdefaultErr();
        boolean ksolveDefault         = ph.getKsolve();
        double kdefaultDefault        = ph.getKdefault();
        double kdefaultErrDefault     = ph.getKdefaultErr();
        boolean debugDefault          = ph.getDebug();
        int verboseDefault            = ph.getVerbose();
		
        // Instantiate a CmdLineParser to read in the arguments pass to PhotomEqSolverRun2b...
        CmdLineParser parser = new CmdLineParser();
        CmdLineParser.Option urlOption             = parser.addStringOption("url");
        CmdLineParser.Option dbNameOption          = parser.addStringOption("dbName");
        CmdLineParser.Option userOption            = parser.addStringOption('u', "user");
        CmdLineParser.Option passwdOption          = parser.addStringOption('p', "passwd");
        CmdLineParser.Option instrumentOption      = parser.addStringOption('i', "instrument");
        CmdLineParser.Option telescopeOption       = parser.addStringOption('t', "telescope");
        CmdLineParser.Option niteOption            = parser.addStringOption('n', "nite");
        CmdLineParser.Option filterOption          = parser.addStringOption('f', "filter");
        CmdLineParser.Option stdColor0Option       = parser.addDoubleOption("stdColor0");
        CmdLineParser.Option ccdidOption           = parser.addIntegerOption('c', "ccdid");
        CmdLineParser.Option magLoOption           = parser.addDoubleOption("magLo");
        CmdLineParser.Option magHiOption           = parser.addDoubleOption("magHi");
        CmdLineParser.Option niterOption           = parser.addIntegerOption("niter");
        CmdLineParser.Option nsigmaOption          = parser.addDoubleOption("nsigma");
        CmdLineParser.Option imageTypeOption       = parser.addStringOption("imageType");
        CmdLineParser.Option imageNameFilterOption = parser.addStringOption("imageNameFilter");
        CmdLineParser.Option runiddescOption       = parser.addStringOption("runiddesc");
        CmdLineParser.Option psmVersionOption      = parser.addStringOption("psmVersion");
        CmdLineParser.Option bsolveOption          = parser.addBooleanOption("bsolve");
        CmdLineParser.Option bdefaultOption        = parser.addDoubleOption("bdefault");
        CmdLineParser.Option bdefaultErrOption     = parser.addDoubleOption("bdefaultErr");
        CmdLineParser.Option ksolveOption          = parser.addBooleanOption("ksolve");
        CmdLineParser.Option kdefaultOption        = parser.addDoubleOption("kdefault");
        CmdLineParser.Option kdefaultErrOption     = parser.addDoubleOption("kdefaultErr");
    	CmdLineParser.Option debugOption           = parser.addBooleanOption('d', "debug");
    	CmdLineParser.Option verboseOption         = parser.addIntegerOption('v', "verbose");
    	CmdLineParser.Option helpOption            = parser.addBooleanOption('h', "help");
    	CmdLineParser.Option paramFileOption       = parser.addStringOption("paramFile");


    	// Process any arguments passed to the main method...
    	try {
    		parser.parse(args);
    	}
    	catch ( CmdLineParser.OptionException e ) {
    		System.err.println(e.getMessage());
    		printUsage(true);
    		System.exit(2);
    	}
    	
    	// If -h or --help was indicated, print help and exit...
    	Boolean help = (Boolean)parser.getOptionValue(helpOption, Boolean.FALSE);
    	if (help) {
    		printUsage(false);
    		System.exit(2);
    	}

    	// If --paramFile option was specified, read the parameter file given...
    	
    	String paramFileName = (String)parser.getOptionValue(paramFileOption);
    	if (paramFileName != null) {

    		System.out.println("paramFile=" + paramFileName);
    		File paramFile = new File(paramFileName);
 
    		if (paramFile.exists() == false || paramFile.canRead() == false) {
    			System.out.println(paramFileName + " either does not exist or cannot be read");
    			System.exit(2);
    		}
    	
    		FileReader fileReader = new FileReader(paramFile);
    		BufferedReader reader = new BufferedReader(fileReader);
    		
    		int iLine = 0;
    		String line = null;
    		
    		//Read in ascii text file containing all the unique star matches in all the overlap regions...
    		System.out.println("Reading contents of parameter file " + paramFileName + ":");
    		while ((line = reader.readLine()) != null) {
    			
    			if (line.length() >= 1 && line.charAt(0) != '#' ) {
     				System.out.println(line);
     				
     				StringTokenizer st = new StringTokenizer(line);
     				int nTokens = st.countTokens();
     				if (nTokens >= 2) {
     					String field1 = st.nextToken();
     					String field2 = st.nextToken();
     					String field = st.toString();
     					if (field1.equals("url")) {
     						urlDefault = field2;
     					} else if (field1.equals("dbName")) {
     						dbNameDefault = field2;
     					} else if (field1.equals("instrument")) {
     						instrumentDefault = field2;
     					} else if (field1.equals("telescope")) {
     						telescopeDefault = field2;
     						System.out.println("telescope = " + field2);
     					} else if (field1.equals("nite")) {
     						niteDefault = field2;
     					} else if (field1.equals("filter")) {
     						filterDefault = field2;
     					} else if (field1.equals("stdColor0")) {
     						stdColor0Default = Double.parseDouble(field2);
     					} else if (field1.equals("ccdid")) {
     						ccdidDefault = Integer.parseInt(field2);
     					} else if (field1.equals("magLo")) {
     						magLoDefault = Double.parseDouble(field2);
     					} else if (field1.equals("magHi")) {
     						magHiDefault = Double.parseDouble(field2);
     					} else if (field1.equals("niterations")) {
     						niterationsDefault = Integer.parseInt(field2);
     					} else if (field1.equals("nsigma")) {
     						nsigmaDefault = Double.parseDouble(field2);
     					} else if (field1.equals("imageType")) {
     						imageTypeDefault = field2;
     					} else if (field1.equals("imageNameFilter")) {
     						imageNameFilterDefault = field2;
     					} else if (field1.equals("runiddesc")) {
     						runiddescDefault = field2;
     					} else if (field1.equals("psmVersion")) {
     						psmVersionDefault = field2;
     					} else if (field1.equals("bsolve")) {
     						bsolveDefault = Boolean.parseBoolean(field2); 
     					} else if (field1.equals("bdefault")) {
     						bdefaultDefault = Double.parseDouble(field2);
     					} else if (field1.equals("bdefaultErr")) {
     						bdefaultErrDefault = Double.parseDouble(field2);
     					} else if (field1.equals("ksolve")) {
     						ksolveDefault = Boolean.parseBoolean(field2);
     					} else if (field1.equals("kdefault")) {
     						kdefaultDefault = Double.parseDouble(field2);
     					} else if (field1.equals("kdefaultErr")) {
     						kdefaultErrDefault = Double.parseDouble(field2);
     					} else if (field1.equals("debug")) {
     						debugDefault = Boolean.parseBoolean(field2);
     					} else if (field1.equals("verbose")) {
     						verboseDefault = Integer.parseInt(field2);
     					}
     					
     				}

     				iLine++;

    			}
    			   			
    		}
    		
    	} else {
    		
    		System.out.println("paramFile=");
    		
    	}

        // Extract values for different parameters; use default values as necessary...
    	String url = (String)parser.getOptionValue(urlOption, urlDefault);
    	String dbName = (String)parser.getOptionValue(dbNameOption, dbNameDefault);
    	String user = (String)parser.getOptionValue(userOption, userDefault);
    	String passwd = (String)parser.getOptionValue(passwdOption, passwdDefault);
    	String instrument = (String)parser.getOptionValue(instrumentOption, instrumentDefault);
    	String telescope = (String)parser.getOptionValue(telescopeOption, telescopeDefault);
    	String nite = (String)parser.getOptionValue(niteOption, niteDefault);
    	String filter = (String)parser.getOptionValue(filterOption, filterDefault);
    	double stdColor0 = ((Double)parser.getOptionValue(stdColor0Option, new Double(stdColor0Default))).doubleValue();
    	int ccdid = ((Integer)parser.getOptionValue(ccdidOption, new Integer(ccdidDefault))).intValue();
    	double magLo = ((Double)parser.getOptionValue(magLoOption, new Double(magLoDefault))).doubleValue();
    	double magHi = ((Double)parser.getOptionValue(magHiOption, new Double(magHiDefault))).doubleValue();
    	int niter = ((Integer)parser.getOptionValue(niterOption, new Integer(niterationsDefault))).intValue();
    	double nsigma = ((Double)parser.getOptionValue(nsigmaOption, new Double(nsigmaDefault))).doubleValue();
    	String imageType = (String)parser.getOptionValue(imageTypeOption, imageTypeDefault);
    	String imageNameFilter = (String)parser.getOptionValue(imageNameFilterOption, imageNameFilterDefault);
    	String runiddesc = (String)parser.getOptionValue(runiddescOption, runiddescDefault);
    	String psmVersion = (String)parser.getOptionValue(psmVersionOption, psmVersionDefault);
    	//Boolean bsolve = (Boolean)parser.getOptionValue(bsolveOption, Boolean.FALSE);
    	//Boolean ksolve = (Boolean)parser.getOptionValue(ksolveOption, Boolean.FALSE);
    	//Boolean debug = (Boolean)parser.getOptionValue(debugOption, Boolean.FALSE);
    	Boolean bsolve = (Boolean)parser.getOptionValue(bsolveOption, bsolveDefault);
    	double bdefault = ((Double)parser.getOptionValue(bdefaultOption, new Double(bdefaultDefault))).doubleValue();
    	double bdefaultErr = ((Double)parser.getOptionValue(bdefaultErrOption, new Double(bdefaultErrDefault))).doubleValue();
    	Boolean ksolve = (Boolean)parser.getOptionValue(ksolveOption, ksolveDefault);
    	double kdefault = ((Double)parser.getOptionValue(kdefaultOption, new Double(kdefaultDefault))).doubleValue();
    	double kdefaultErr = ((Double)parser.getOptionValue(kdefaultErrOption, new Double(kdefaultErrDefault))).doubleValue();
    	Boolean debug = (Boolean)parser.getOptionValue(debugOption, debugDefault);
    	int verbose = ((Integer)parser.getOptionValue(verboseOption, new Integer(verboseDefault))).intValue();
    	
    	
    	// Grab any other arguments that were not passed as part of an option/value pair...
    	// (there shouldn't be any such arguments, but check anyway...)
    	String[] otherArgs = parser.getRemainingArgs();
    	
    	
    	// Set the instance variables for the PhotomEqSolver2b object ph using the values 
    	// determined above...
    	System.out.println("\n\nSetting the values of the paramters:");

    	ph.setUrl(url);   
    	System.out.println("url="+ph.getUrl());
    	
    	ph.setDbName(dbName);   
    	System.out.println("dbName="+ph.getDbName());
    	
    	ph.setUser(user);   
    	System.out.println(user);
    	
    	ph.setPasswd(passwd);   
    	System.out.println(passwd);
    	
    	ph.setInstrument(instrument);   
    	System.out.println("instrument="+ph.getInstrument());
    	
    	ph.setTelescope(telescope);   
    	System.out.println("telescope="+ph.getTelescope());
    	
    	ph.setNite(nite);
    	System.out.println("nite="+ph.getNite());
    	
    	ph.setFilter(filter);
    	System.out.println("filter="+ph.getFilter());
    	
    	ph.setStdColor0(stdColor0);
    	System.out.println("stdColor0="+ph.getStdColor0());
    	
    	ph.setCcdid(ccdid);
    	System.out.println("ccdid="+ph.getCcdid());
    	
    	ph.setMagLo(magLo);
    	System.out.println("magLo="+ph.getMagLo());
    	
    	ph.setMagHi(magHi);
    	System.out.println("magHi="+ph.getMagHi());
    	
    	ph.setNiterations(niter);
    	System.out.println("niterations="+ph.getNiterations());
    	
    	ph.setNsigma(nsigma);   
    	System.out.println("nsigma="+ph.getNsigma());
    	
    	ph.setImageType(imageType);   
    	System.out.println("imageType="+ph.getImageType());
    	
    	ph.setImageNameFilter(imageNameFilter);   
    	System.out.println("imageNameFilter="+ph.getImageNameFilter());
    	
    	ph.setRuniddesc(runiddesc);   
    	System.out.println("runiddesc="+ph.getRuniddesc());
    	
    	ph.setPsmVersion(psmVersion);   
    	System.out.println("psmVersion="+ph.getPsmVersion());
       	
    	ph.setBsolve(bsolve);
    	System.out.println("bsolve="+ph.getBsolve());
       	
    	ph.setBdefault(bdefault);
    	System.out.println("bdefault="+ph.getBdefault());
       	
    	ph.setBdefaultErr(bdefaultErr);
    	System.out.println("bdefaultErr="+ph.getBdefaultErr());
    	
    	ph.setKsolve(ksolve);
    	System.out.println("ksolve="+ph.getKsolve());
       	
    	ph.setKdefault(kdefault);
    	System.out.println("kdefault="+ph.getKdefault());
       	
    	ph.setKdefaultErr(kdefaultErr);
    	System.out.println("kdefaultErr="+ph.getKdefaultErr());
    	
    	ph.setDebug(debug);
    	System.out.println("debug="+ph.getDebug());
    	
    	ph.setVerbose(verbose);   
    	System.out.println("verbose="+ph.getVerbose());
    	
    	
    	// These are instance variables that should not change, so
    	// there are no command line options available to them.
    	// Nonetheless, we set them here...
    	ph.setSqlDriver("oracle.jdbc.driver.OracleDriver");
    	ph.setStdTable("standard_stars");
    	ph.setObsTable("OBJECTS");
    	ph.setFilesTable("FILES");
    	ph.setFitTable("psmfit");
    	
    	Date date = new Date();
    	ph.setDate(date);
    	
    	if (true) {
    		try {
    			ph.solve();
    		} catch (ClassNotFoundException e) {
    			e.printStackTrace();
    		} catch (SQLException e) {
    			e.printStackTrace();
    		} catch (Exception e) {
    			e.printStackTrace();
    		}
    	}
    	
    }
    
}
