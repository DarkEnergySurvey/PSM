package gov.fnal.eag.dtucker.desPhotoStds.test;

import java.io.FileNotFoundException;
import java.io.IOException;

import gov.fnal.eag.dtucker.desPhotoStds.DESDMFile;
import junit.framework.TestCase;

public class DESDMFileTest extends TestCase {

	/*
	 * Test method for 'gov.fnal.eag.dtucker.util.DESDMFile.readDESDMFile()'
	 */
	public void testReadDESDMFile() throws IOException, ClassNotFoundException, FileNotFoundException {
		System.out.println("testReadDESDMFile...");
		DESDMFile desdmFile = new DESDMFile();
		desdmFile.readDESDMFile();
		System.out.println(
				"DB_USER:  \t" + desdmFile.getDbUser() + "\n" + 
				"DB_PASSWD:\t" + desdmFile.getDbPasswd() + "\n" + 
				"DB_SERVER:\t" + desdmFile.getDbServer() + "\n" + 
				"DB_NAME:  \t" + desdmFile.getDbName() + "\n");
	}
	
	public void testToAsterisks() {

		System.out.println("testToAsterisks...");
		DESDMFile desdmFile = new DESDMFile();
		
		String testString;
		String testStringAsterisks;

		testString = "a";
		testStringAsterisks = desdmFile.toAsterisks(testString);
		System.out.println(testString + ":\t\t\t" + testStringAsterisks);
		assertEquals(testStringAsterisks, "*");
		
		testString = "abc";
		testStringAsterisks = desdmFile.toAsterisks(testString);
		System.out.println(testString + ":\t\t\t" + testStringAsterisks);
		assertEquals(testStringAsterisks, "***");
		
		testString = "abc123";
		testStringAsterisks = desdmFile.toAsterisks(testString);
		System.out.println(testString + ":\t\t\t" + testStringAsterisks);
		assertEquals(testStringAsterisks, "******");
	
		testString = "abc 123";
		testStringAsterisks = desdmFile.toAsterisks(testString);
		System.out.println(testString + ":\t\t" + testStringAsterisks);
		assertEquals(testStringAsterisks, "*******");
	
		testString = "abc 123 @ xyz : ijk";
		testStringAsterisks = desdmFile.toAsterisks(testString);
		System.out.println(testString + ":\t" + testStringAsterisks);
		assertEquals(testStringAsterisks, "*******************");	

		System.out.println("");
		
	}

}
