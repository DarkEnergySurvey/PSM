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
public class ReadAnASCIIFile {

    public static void main(String[] args)  {
        int nargs = args.length;
        String fileName;
        if (nargs > 0) {
            fileName = args[0]; 
        } else {
            fileName = "WriteAFile.txt";
        }
        try {
            File myFile = new File(fileName);
            if (myFile.exists() && myFile.canRead()) {
                FileReader fileReader = new FileReader(myFile);
                BufferedReader reader = new BufferedReader(fileReader);
                String line = null;
                while ((line = reader.readLine()) != null) {
                    System.out.println(line);
                }
                reader.close();
            } else {
                System.out.println(fileName + " either does not exist or cannot be read");
            }
        } catch(Exception ex) {
                ex.printStackTrace();
        }
    }
}
