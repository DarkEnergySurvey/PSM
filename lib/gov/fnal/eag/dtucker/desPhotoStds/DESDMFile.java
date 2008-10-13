package gov.fnal.eag.dtucker.desPhotoStds;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

public class DESDMFile {

	private static String desdmFileName = ".desdm";
	private String dbUser = "";
	private String dbPasswd = "";
	private String dbServer = "";
	private String dbName = "";
	
	public void readDESDMFile() throws IOException, ClassNotFoundException, FileNotFoundException {
		
		File desdmFile = new File(System.getProperty("user.home"),desdmFileName);
		String desdmFullFileName = desdmFile.getCanonicalPath();
		
		if (!desdmFile.exists()) {
			System.err.println(desdmFullFileName + " does not exist!");
			throw new FileNotFoundException();
		}
		
		if (!desdmFile.canRead()) {
			System.err.println(desdmFullFileName + " can not be read!");
			throw new IOException();
		}
		
		FileReader fileReader = new FileReader(desdmFile);
		BufferedReader reader = new BufferedReader(fileReader);
		String line = null;
		
		while ((line = reader.readLine()) != null) {
			
			if (line.length() == 0) {
				continue;
			}
			if (line.charAt(0) == '#') {
				continue;
			}
			
			StringTokenizer st = new StringTokenizer(line);
			int nTokens = st.countTokens();
			if (nTokens != 2) {
				continue;
			}
			
			String param = st.nextToken();
			String value = st.nextToken();
			
			if (param.equals("DB_USER")) {
				this.setDbUser(value);
			} else if (param.equals("DB_PASSWD")) {
				this.setDbPasswd(value);
			} else if (param.equals("DB_SERVER")) {
				this.setDbServer(value);
			} else if (param.equals("DB_NAME")) {
				this.setDbName(value);
			}
			
		}
		
		reader.close();
		
	}
	
	/**
	 * A method to hide the values of dbUser and dbPasswd behind a string of asterisks.
	 * @param string
	 * @return asterisks
	 */
	public String toAsterisks(String string) {
		int length = string.length();
		String asterisks = "";
		for (int i=0;i<length;i++) {
			asterisks = asterisks + "*";
		}
		return asterisks;
	}

	public String getDbName() {
		return dbName;
	}

	public void setDbName(String dbName) {
		this.dbName = dbName;
	}

	public String getDbPasswd() {
		return dbPasswd;
	}

	public void setDbPasswd(String dbPasswd) {
		this.dbPasswd = dbPasswd;
	}

	public String getDbServer() {
		return dbServer;
	}

	public void setDbServer(String dbServer) {
		this.dbServer = dbServer;
	}

	public String getDbUser() {
		return dbUser;
	}

	public void setDbUser(String dbUser) {
		this.dbUser = dbUser;
	}		

}





