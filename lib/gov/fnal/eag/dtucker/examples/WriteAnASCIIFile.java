/*
 * Created on Aug 11, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
package gov.fnal.eag.dtucker.examples;
import java.io.*;


/**
 * @author dtucker
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
public class WriteAnASCIIFile {

    public static void main(String[] args) {
        try {
            FileWriter writer = new FileWriter("WriteAFile.txt");
            writer.write("how about this?\n");
            writer.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
}
