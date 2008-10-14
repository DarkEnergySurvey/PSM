/**
 * PhotomEqSolverRunDC4.java   The DES Collaboration, 14 September 2008
 */

package gov.fnal.eag.dtucker.desPhotoStds;

import gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverDC4;

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
public class PhotomEqSolverRunDC4 {
	
	private static void printUsage(boolean error) {

		PhotomEqSolverDC4 ph = new PhotomEqSolverDC4();
        
		// Grab default values from PhotomEqSolverDC4...
		String urlDefault             = ph.getUrl();
        String dbNameDefault          = ph.getDbName();
        String userDefault            = ph.getUser();
        String passwdDefault          = ph.getPasswd();
        String projectDefault         = ph.getProject();
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
        String runDefault             = ph.getRun();
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
       			"\nUsage:  java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRunDC4 [OPTIONS] \n\n" +
				"   OPTIONS: \n" + 
				"   --paramFile VALUE             name of (optional) parameter file \n" + 
				"                                 values for parameters in the paramFile will override values hardwired into the PSMDC4 code \n" +
				"                                 (or in a .desdm file), but values passed as optional parameters in the command line will \n" +
				"                                 override the values in the paramFile (or in a .desdm file) \n" +
				"   --url VALUE                   URL of database                       [default: " + urlDefault + "] \n" + 
				"   --dbName VALUE                name of database                      [default: " + dbNameDefault + "] \n" + 
				"   -u VALUE, --user VALUE        username for db access                [default: " + userDefault + "] \n" + 
				"   -p VALUE, --passwd VALUE      password for db access                [default: " + passwdDefault + "] \n" + 
				"   -P VALUE, --project VALUE     name of project in db                 [default: " + projectDefault + "] \n" +
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
				"   --run VALUE                   run for std star fields               [default: " + runDefault + "] \n" + 
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
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRunDC4 --url jdbc:oracle:thin:@charon.ncsa.uiuc.edu:1521: --dbName des  -u myUserName -p myPassword -P BCS -n 20061223 -f g --ccdid 0 --magLo 15.0 --magHi 18.0 --niter 3 --nsigma 2.5 --imageType remap --imageNameFilter % --run 20080324000000_20061223 --psmVersion v_DC4 --debug -v 2 --bsolve --ksolve \n\n" +
				"   Example 1a (with .desdm file): \n" +
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRunDC4 -P BCS -n 20061223 -f g --ccdid 0 --magLo 15.0 --magHi 18.0 --niter 3 --nsigma 2.5 --imageType remap --imageNameFilter % --run 20080324000000_20061223 --psmVersion v_DC4 --debug -v 2 --bsolve --ksolve \n\n" +
				"   Example 2: \n" +
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRunDC4 --user myUserName -p myPassword -P BCS -n 20061223 -f g --debug -v 2 --ksolve \n\n" + 
				"   Example 2a (with .desdm file): \n" +
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRunDC4 -P BCS -n 20061223 -f g --debug -v 2 --ksolve \n\n" + 
				"   Example 3: \n" +
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRunDC4 --user myUserName -p myPassword -P BCS --paramFile /home/myname/psmParam.par \n\n" + 
				"   Example 3a (with .desdm file): \n" +
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRunDC4 -P BCS --paramFile /home/myname/psmParam.par \n\n" + 
				"   Example 4: \n" +
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRunDC4 --user myUserName -p myPassword -P BCS -n 20061223 -f g --debug -v 2 --ksolve  --paramFile /home/myname/psmParam.par \n\n" + 
				"   Example 4a (with .desdm file): \n" +
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRunDC4 -P BCS -n 20061223 -f g --debug -v 2 --ksolve  --paramFile /home/myname/psmParam.par \n\n" + 
				"   Example 5 (no .desdm file required): \n" + 
				"      java gov.fnal.eag.dtucker.desPhotoStds.PhotomEqSolverRunDC4 --help \n\n";


				
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

    	int localVerbose = 0;
        for (int i=0; i< args.length; i++) {
        	if (args[i].equals("-v") || args[i].equals("--verbose")) {
        		if (i+1 < args.length) {
        			localVerbose = Integer.parseInt(args[i+1]);
        		}
        	}
        }

    	boolean ignoreParamFileBTermInfo  = false;
        for (int i=0; i< args.length; i++) {
        	if (args[i].equals("--bdefault") || args[i].equals("--bdefaultErr")) {
        		ignoreParamFileBTermInfo = true;
        		break;
        	}
        }
        //System.out.println("ignoreParamFileBTermInfo="+ignoreParamFileBTermInfo);
    	    	
        if (localVerbose > 0) {
        	System.out.println("PhotomEqSolverRunDC4 \n");
        }
 
		if (localVerbose > 0) {
			System.out.println("localVerbose = " + localVerbose  + "\n");
		}
        
        if (localVerbose > 0) {
        	System.out.print("arglist:  ");
        	for (int i=0; i< args.length; i++) {
        		// protect against displaying user name or password in the clear...
        		if (	( args[i].equals("-u") || args[i].equals("--user") ||
        				  args[i].equals("-p") || args[i].equals("--passwd") )  && 
        				( i+1 < args.length )    ) {
        			System.out.print(args[i] + " ");
        			System.out.print(new DESDMFile().toAsterisks(args[i+1]) + " ");
        			i++;
        		} else {
        			System.out.print(args[i] + " ");
        		}
        	}
        	System.out.print("\n\n");
        }
        	
        // Instantiate an instance of the PhotomEqSolverDC4 class...
        PhotomEqSolverDC4 ph = new PhotomEqSolverDC4();

        // Grab default values for parameters from PhotomEqSolverDC4...
		String urlDefault             = ph.getUrl();
        String dbNameDefault          = ph.getDbName();
        String userDefault            = ph.getUser();
        String passwdDefault          = ph.getPasswd();
        String projectDefault         = ph.getProject();
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
        String runDefault             = ph.getRun();
        String psmVersionDefault      = ph.getPsmVersion();
        boolean bsolveDefault         = ph.getBsolve();
        double bdefaultDefault        = ph.getBdefault();
        double bdefaultErrDefault     = ph.getBdefaultErr();
        boolean ksolveDefault         = ph.getKsolve();
        double kdefaultDefault        = ph.getKdefault();
        double kdefaultErrDefault     = ph.getKdefaultErr();
        boolean debugDefault          = ph.getDebug();
        int verboseDefault            = ph.getVerbose();

        // Instantiate an instance of the ColorTermCoeffs class...
        ColorTermCoeffs colorTermCoeffs = new ColorTermCoeffs();
		
        
		// Grab default values for user, passwd, url, and dbName from .desdm (if it exists).
        // These values take precedence over those hardwired in the PhotomEqSolverDC4 code, 
        // but can be overwritten by values in the paramFile (url, dbName) or in the 
        // argument list (user, passwd, url, dbName)...
		DESDMFile desdmFile = new DESDMFile();
		try {
			if (localVerbose > 0) {
				System.out.println("Reading .desdm file... \n");
			}
			desdmFile.readDESDMFile();
			userDefault = desdmFile.getDbUser();
			passwdDefault = desdmFile.getDbPasswd();
			dbNameDefault = desdmFile.getDbName();
			urlDefault = "jdbc:oracle:thin:@" + desdmFile.getDbServer() + ":1521:";
		} catch (Exception e) {
			if (localVerbose > 0) {
				System.err.println("Problem accessing or reading .desdm file.  Skipping... \n");
			}
		}

		// Instantiate a CmdLineParser to read in the arguments pass to PhotomEqSolverRunDC4...
        CmdLineParser parser = new CmdLineParser();
        CmdLineParser.Option urlOption             = parser.addStringOption("url");
        CmdLineParser.Option dbNameOption          = parser.addStringOption("dbName");
        CmdLineParser.Option userOption            = parser.addStringOption('u', "user");
        CmdLineParser.Option passwdOption          = parser.addStringOption('p', "passwd");
        CmdLineParser.Option projectOption         = parser.addStringOption('P', "project");
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
        CmdLineParser.Option runOption             = parser.addStringOption("run");
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

    		if (localVerbose > 1) {
    			System.out.println("paramFile=" + paramFileName + "\n");
    		}
    		File paramFile = new File(paramFileName);
 
    		if (paramFile.exists() == false || paramFile.canRead() == false) {
    			if (localVerbose > 0) {
    				System.out.println(paramFileName + " either does not exist or cannot be read \n");
    			}
    			System.exit(2);
    		}
    	
    		FileReader fileReader = new FileReader(paramFile);
    		BufferedReader reader = new BufferedReader(fileReader);
    		
    		int iLine = 0;
    		String line = null;
    		
    		//Read in parameter file...
    		if (localVerbose > 0) {
    			System.out.println("Reading contents of parameter file " + paramFileName + "...");
    		}
    		while ((line = reader.readLine()) != null) {
    			
    			if (line.length() >= 1 && line.charAt(0) != '#' ) {
    				
    				if (localVerbose > 1) {
    					System.out.println(line);
    				}
    				
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
     					} else if (field1.equals("project")) {
     						projectDefault = field2;
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
     					} else if (field1.equals("run")) {
     						runDefault = field2;
     					} else if (field1.equals("psmVersion")) {
     						psmVersionDefault = field2;
     					} else if (field1.equals("bsolve")) {
     						bsolveDefault = Boolean.parseBoolean(field2); 
     					//} else if (field1.equals("bdefault")) {
     					//	bdefaultDefault = Double.parseDouble(field2);
     					//} else if (field1.equals("bdefaultErr")) {
     					//	bdefaultErrDefault = Double.parseDouble(field2);
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
     					} else if (field1.equals("bccdidArray")) {
     						
     						if (ignoreParamFileBTermInfo == false) {
     							ArrayList tempArrayList = new ArrayList();
     							//System.out.print("bccdidArray: ");
     							//System.out.print(field2 + " ");
     							tempArrayList.add(field2);
     							for (int iToken = 2; iToken < nTokens; iToken++) {
     								String fieldi = st.nextToken();
     								if (fieldi.charAt(0) == '#') {break;}
     								//System.out.print(fieldi + " ");
     								tempArrayList.add(fieldi);
     							}
     							//System.out.println("");
     							colorTermCoeffs.setBccdidArrayList(tempArrayList);
     						
     							if (false) {
     								System.out.print("bccdidArray: ");
     								for (int jj=0;jj<colorTermCoeffs.getBccdidArrayList().size();jj++) {
     									System.out.print(colorTermCoeffs.getBccdidArrayList().get(jj) + " ");
     								}
     								System.out.println("");
     							}
     							
     						}
     						
     					} else if (field1.equals("bdefaultArray")) {
     						
     						if (ignoreParamFileBTermInfo == false) {
     							ArrayList tempArrayList = new ArrayList();
     							//System.out.print("bdefaultArray: ");
     							//System.out.print(field2 + " ");
     							tempArrayList.add(field2);
     							for (int iToken = 2; iToken < nTokens; iToken++) {
     								String fieldi = st.nextToken();
     								if (fieldi.charAt(0) == '#') {break;}
     								//System.out.print(fieldi + " ");
     								tempArrayList.add(fieldi);
     							}
     							//System.out.println("");
     							colorTermCoeffs.setBdefaultArrayList(tempArrayList);
     							
     							if (false) {
     								System.out.print("bdefaultArray: ");
     								for (int jj=0;jj<colorTermCoeffs.getBdefaultArrayList().size();jj++) {
     									System.out.print(colorTermCoeffs.getBdefaultArrayList().get(jj) + " ");
     								}
     								System.out.println("");
     							}
     							
     						}
     						
     					} else if (field1.equals("bdefaultErrArray")) {
     						
     						if (ignoreParamFileBTermInfo == false) {
     							ArrayList tempArrayList = new ArrayList();
     							//System.out.print("bdefaultErrArray: ");
     							//System.out.print(field2 + " ");
     							tempArrayList.add(field2);
     							for (int iToken = 2; iToken < nTokens; iToken++) {
     								String fieldi = st.nextToken();
     								if (fieldi.charAt(0) == '#') {break;}
     								//System.out.print(fieldi + " ");
     								tempArrayList.add(fieldi);
     							}
     							//System.out.println("");
     							colorTermCoeffs.setBdefaultErrArrayList(tempArrayList);
     							
     							if (false) {
     								System.out.print("bdefaultErrArray: ");
     								for (int jj=0;jj<colorTermCoeffs.getBdefaultErrArrayList().size();jj++) {
     									System.out.print(colorTermCoeffs.getBdefaultErrArrayList().get(jj) + " ");
     								}
     								System.out.println("");
     							}
     							
     						}     						
     						
     					}
     					
     				}

     				iLine++;

    			}
    			   			
    		}
    		
    	} else {
    		
    		if (localVerbose > 0) {
    			System.out.println("paramFile=\n");
    		}
    		
    	}

        // Extract values for different parameters; use default values as necessary...
    	String url = (String)parser.getOptionValue(urlOption, urlDefault);
    	String dbName = (String)parser.getOptionValue(dbNameOption, dbNameDefault);
    	String user = (String)parser.getOptionValue(userOption, userDefault);
    	String passwd = (String)parser.getOptionValue(passwdOption, passwdDefault);
    	String project = (String)parser.getOptionValue(projectOption, projectDefault);
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
    	String run = (String)parser.getOptionValue(runOption, runDefault);
    	String psmVersion = (String)parser.getOptionValue(psmVersionOption, psmVersionDefault);
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
    	
    	
    	// if the ArrayLists in colorTermCoeffs are empty, set to default values...
    	if (colorTermCoeffs.getBccdidArrayList().size() == 0 ||
    			colorTermCoeffs.getBdefaultArrayList().size() == 0 ||
    			colorTermCoeffs.getBdefaultErrArrayList().size() == 0 ) {

    		ArrayList tempBccdidArrayList = new ArrayList();
    		tempBccdidArrayList.add(new Integer (ccdid));
    		colorTermCoeffs.setBccdidArrayList(tempBccdidArrayList);

    		ArrayList tempBdefaultArrayList = new ArrayList();
    		tempBdefaultArrayList.add(new Double (bdefault));
    		colorTermCoeffs.setBdefaultArrayList(tempBdefaultArrayList);
    	
    		ArrayList tempBdefaultErrArrayList = new ArrayList();
    		tempBdefaultErrArrayList.add(new Double (bdefaultErr));
    		colorTermCoeffs.setBdefaultErrArrayList(tempBdefaultErrArrayList);

    	}
    	
    	
    	// Set the instance variables for the PhotomEqSolverDC4 object ph using the values 
    	// determined above...
    	if (localVerbose > 0) {
    		System.out.println("\n\nSetting the values of the parameters:");
    	}
    	
    	ph.setUrl(url);   
    	if (localVerbose > 0) {System.out.println("url="+ph.getUrl());}
    	
    	ph.setDbName(dbName);   
    	if (localVerbose > 0) {System.out.println("dbName="+ph.getDbName());}
    	
    	ph.setUser(user);   
    	if (localVerbose > 0) {System.out.println("user="+desdmFile.toAsterisks(ph.getUser()));}
    	
    	ph.setPasswd(passwd);   
    	if (localVerbose > 0) {System.out.println("passwd="+desdmFile.toAsterisks(ph.getPasswd()));}
    	
    	ph.setProject(project);   
    	if (localVerbose > 0) {System.out.println("project="+ph.getProject());}
    	
    	ph.setNite(nite);
    	if (localVerbose > 0) {System.out.println("nite="+ph.getNite());}
    	
    	ph.setFilter(filter);
    	if (localVerbose > 0) {System.out.println("filter="+ph.getFilter());}
    	
    	ph.setStdColor0(stdColor0);
    	if (localVerbose > 0) {System.out.println("stdColor0="+ph.getStdColor0());}
    	
    	ph.setCcdid(ccdid);
    	if (localVerbose > 0) {System.out.println("ccdid="+ph.getCcdid());}
    	
    	ph.setMagLo(magLo);
    	if (localVerbose > 0) {System.out.println("magLo="+ph.getMagLo());}
    	
    	ph.setMagHi(magHi);
    	if (localVerbose > 0) {System.out.println("magHi="+ph.getMagHi());}
    	
    	ph.setNiterations(niter);
    	if (localVerbose > 0) {System.out.println("niterations="+ph.getNiterations());}
    	
    	ph.setNsigma(nsigma);   
    	if (localVerbose > 0) {System.out.println("nsigma="+ph.getNsigma());}
    	
    	ph.setImageType(imageType);   
    	if (localVerbose > 0) {System.out.println("imageType="+ph.getImageType());}
    	
    	ph.setImageNameFilter(imageNameFilter);   
    	if (localVerbose > 0) {System.out.println("imageNameFilter="+ph.getImageNameFilter());}
    	
    	ph.setRun(run);   
    	if (localVerbose > 0) {System.out.println("run="+ph.getRun());}
    	
    	ph.setPsmVersion(psmVersion);   
    	if (localVerbose > 0) {System.out.println("psmVersion="+ph.getPsmVersion());}
       	
    	ph.setBsolve(bsolve);
    	if (localVerbose > 0) {System.out.println("bsolve="+ph.getBsolve());}

    	ph.setBdefault(bdefault);
    	if (localVerbose > 0) {System.out.println("bdefault="+ph.getBdefault());}

		ph.setBdefaultErr(bdefaultErr);
		if (localVerbose > 0) {System.out.println("bdefaultErr="+ph.getBdefaultErr());}

    	ph.setColorTermCoeffs(colorTermCoeffs);
    	
    	if (localVerbose > 0) {
    		System.out.print("bccdidArray=");
    		for (int jj=0;jj<ph.getColorTermCoeffs().getBccdidArrayList().size();jj++) {
    			System.out.print(ph.getColorTermCoeffs().getBccdidArrayList().get(jj) + " ");
    		}
    		System.out.println("");    		
    	}
    	
    	if (localVerbose > 0) {
    		System.out.print("bdefaultArray=");
    		for (int jj=0;jj<ph.getColorTermCoeffs().getBdefaultArrayList().size();jj++) {
    			System.out.print(ph.getColorTermCoeffs().getBdefaultArrayList().get(jj) + " ");
    		}
    		System.out.println("");    		
    	}
       	
    	if (localVerbose > 0) {
    		System.out.print("bdefaultErrArray=");
    		for (int jj=0;jj<ph.getColorTermCoeffs().getBdefaultErrArrayList().size();jj++) {
    			System.out.print(ph.getColorTermCoeffs().getBdefaultErrArrayList().get(jj) + " ");
    		}
    		System.out.println("");
    	}
    	
    	ph.setKsolve(ksolve);
    	if (localVerbose > 0) {System.out.println("ksolve="+ph.getKsolve());}
       	
    	ph.setKdefault(kdefault);
    	if (localVerbose > 0) {System.out.println("kdefault="+ph.getKdefault());}
       	
    	ph.setKdefaultErr(kdefaultErr);
    	if (localVerbose > 0) {System.out.println("kdefaultErr="+ph.getKdefaultErr());}
    	
    	ph.setDebug(debug);
    	if (localVerbose > 0) {System.out.println("debug="+ph.getDebug());}
    	
    	ph.setVerbose(verbose);   
    	if (localVerbose > 0) {System.out.println("verbose="+ph.getVerbose());}
    	
    	
    	// These are instance variables that should not change, so
    	// there are no command line options available to them.
    	// Nonetheless, we set them here...
    	ph.setSqlDriver("oracle.jdbc.driver.OracleDriver");
    	ph.setStdTable("standard_stars");
    	ph.setObsTable("OBJECTS");
    	ph.setImageTable("IMAGE");
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
