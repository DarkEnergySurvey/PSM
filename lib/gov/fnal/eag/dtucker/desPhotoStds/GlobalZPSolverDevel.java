package gov.fnal.eag.dtucker.desPhotoStds;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.Arrays;
import java.util.TreeMap;
import java.util.TreeSet;

import nom.tam.fits.Fits;
import nom.tam.util.BufferedFile;

import cern.colt.list.DoubleArrayList;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;


/**
 * This class defines methods for solving the frame-to-frame zeropoint offsets
 * for a set of overlapping frames
 * @author dtucker
 *
 */
public class GlobalZPSolverDevel {
	
	//Instance variables dealing with this run of the Global Calibrations Module
	private Date date = new Date();
	
	//Instance variables dealing with the observed data to be calibrated
	private String starMatchFileName  = null;
	private String outputFileName = null;
	
	//General instance variables...
	private String[] filterList = {"u", "g", "r", "i", "z"};
	private int verbose = 0;
	
	public void solve () throws Exception {
	
		System.out.println("\n\nGlobalZPSolver");
		System.out.println("Start Time: \t" + (new Date()).toString());

		if (verbose > 0) {
			System.out.println("");
			System.out.println("The beginning...");
			System.out.println("");
		}
		
		
		File starMatchFile = new File(starMatchFileName);
		if (starMatchFile.exists() == false || starMatchFile.canRead() == false) {
			System.out.println(starMatchFileName + " either does not exist or cannot be read");
			return;
		}
		
		FileReader fileReader = new FileReader(starMatchFile);
		BufferedReader reader = new BufferedReader(fileReader);
		
		ArrayList regionidList = new ArrayList();
		ArrayList matchedStarList = new ArrayList();
		ArrayList regionRawOffsetsList = new ArrayList();
		int regionRawOffsetsListSize = 0;
				
		Set regionPairSet = new TreeSet();		
		Map regionPairMap = new TreeMap();
		Map regionQualityMap = new TreeMap();
				
		int iLine = 0;
		String line = null;
		
		//Read in ascii text file containing all the unique star matches in all the overlap regions...
		while ((line = reader.readLine()) != null) {
			
			if (line.length() == 0) {
				System.out.println("line \'" + line + "\' has zero length...  skipping...");
				continue;
			}
			
			if (line.charAt(0) == '#') {
				System.out.println("line \'" + line + "\' commented out...  skipping...");
				continue;
			}
			
			StringTokenizer st = new StringTokenizer(line);
			int nTokens = st.countTokens();
			if (nTokens != 19) {
				System.out.println("line \'" + line + "\' does not contain 19 fields...  skipping...");
				continue;
			}
			
			MatchedStar matchedStar = new MatchedStar();
			
			int regionid1 = Integer.parseInt(st.nextToken());
			matchedStar.setRegionid1(regionid1);
			matchedStar.setRegionRaCenDeg1(Double.parseDouble(st.nextToken()));
			matchedStar.setRegionDecCenDeg1(Double.parseDouble(st.nextToken()));
			int regionQuality1 = Integer.parseInt(st.nextToken());
			matchedStar.setRegionQuality1(regionQuality1);
			matchedStar.setStarid1(Integer.parseInt(st.nextToken()));
			matchedStar.setRaDeg1(Double.parseDouble(st.nextToken()));
			matchedStar.setDecDeg1(Double.parseDouble(st.nextToken()));
			double mag1 = Double.parseDouble(st.nextToken());
			matchedStar.setMag1(mag1);
			matchedStar.setMagErr1(Double.parseDouble(st.nextToken()));

			int regionid2 = Integer.parseInt(st.nextToken());
			matchedStar.setRegionid2(regionid2);
			matchedStar.setRegionRaCenDeg2(Double.parseDouble(st.nextToken()));
			matchedStar.setRegionDecCenDeg2(Double.parseDouble(st.nextToken()));
			int regionQuality2 = Integer.parseInt(st.nextToken());
			matchedStar.setRegionQuality2(regionQuality2);
			matchedStar.setStarid2(Integer.parseInt(st.nextToken()));
			matchedStar.setRaDeg2(Double.parseDouble(st.nextToken()));
			matchedStar.setDecDeg2(Double.parseDouble(st.nextToken()));
			double mag2 = Double.parseDouble(st.nextToken());
			matchedStar.setMag2(mag2);
			matchedStar.setMagErr2(Double.parseDouble(st.nextToken()));

			matchedStar.setSepArcsec12(Double.parseDouble(st.nextToken()));
			
			matchedStarList.add(matchedStar);
			
			double dmag = mag1-mag2;

			//If we haven't encountered this region before, add its id to the regionidList...
			if (regionidList.contains(regionid1) == false) {
				regionidList.add(new Integer(regionid1));
				//if (regionidList.size() == 1) {
					//regionQualityMap.put(regionid1, (Integer) 1);
				//} else {
					//regionQualityMap.put(regionid1, (Integer) 0);
				//}
				if (regionQuality1 == 1) {
					if (Math.random() < 0.90) {	
						regionQuality1 = 0;
						System.out.println(regionid1 + "\t" + regionQuality1);
					}
				}
				regionQualityMap.put(regionid1,(Integer) regionQuality1);
			}
			if (regionidList.contains(regionid2) == false) {
				regionidList.add(new Integer(regionid2));
				//regionQualityMap.put(regionid2, (Integer) 0);
				if (regionQuality2 == 1) {
					if (Math.random() < 0.90) {	
						regionQuality2 = 0;
						System.out.println(regionid2 + "\t" + regionQuality2);
					}
				}
				regionQualityMap.put(regionid2,(Integer) regionQuality2);
			}

			String regionPairName = regionid1 + " " + regionid2;
			regionPairSet.add(regionPairName);			
			
			if (regionPairMap.containsKey(regionPairName) == false) {
				regionPairMap.put(regionPairName, (Integer) regionRawOffsetsListSize);
				RegionRawOffsets rro = new RegionRawOffsets();
				rro.setRegionid1(regionid1);
				rro.setRegionid2(regionid2);
				DoubleArrayList dmagList = new DoubleArrayList();
				dmagList.add(new Double(dmag));
				rro.setDmagList(dmagList);
				regionRawOffsetsList.add(rro);
				regionRawOffsetsListSize++;	
			} else {
				int kk = (Integer) regionPairMap.get(regionPairName);
				RegionRawOffsets rro = (RegionRawOffsets) regionRawOffsetsList.get(kk);
				DoubleArrayList dmagList = rro.getDmagList();
				dmagList.add(new Double(dmag));
				rro.setDmagList(dmagList);		
			}
			
			if (verbose > 0) {
				if (    (iLine <     100) || 
						(iLine <    1000 && iLine%100  == 0) || 
						(iLine <   10000 && iLine%1000 == 0) || 
						(iLine%10000 == 0) 
						) {
					System.out.println("line " + iLine + ":  " + line);
				}
			}
			
			iLine++;
			
		}

		reader.close();
		int matchedStarListSize = matchedStarList.size();

		if (verbose > 0) {
			System.out.println("There are " + matchedStarListSize + " entries in matchedStarList.");
		}
				
		if (verbose > 0) {
			System.out.println("There are " + regionRawOffsetsListSize + " entries in regionRawOffsetsList.");
		}
		
		if (verbose > 1) {
			System.out.println("TreeMap regionPairMap contains " + regionPairMap.size() + " entries");
			System.out.println("output of TreeMap map:");
			Iterator iter = regionPairSet.iterator();
			while (iter.hasNext()) {
				String regionPairName = (String) iter.next();
				System.out.println(regionPairName +  " " + regionPairMap.get(regionPairName));
			}
		}

		//Sort regionidList...
		Collections.sort(regionidList);	
		int regionidListSize = regionidList.size();
		if (verbose > 1) {
			for (int ii=0; ii<regionidListSize; ii++) {
				System.out.println(ii + " " + Integer.parseInt(regionidList.get(ii).toString()));
			}
		}

				
		//Loop through all the regionRawOffsets objects on the regionRawOffsetsList, 
		// and calculate the median, mean, std error of the dmagList for each region1-region2
		// overlap region...
		for (int kk=0; kk<regionRawOffsetsListSize; kk++) {
			RegionRawOffsets rro = (RegionRawOffsets) regionRawOffsetsList.get(kk);
			int regionid1 = rro.getRegionid1();
			int regionid2 = rro.getRegionid2();
			DoubleArrayList dmagList = rro.getDmagList();
			if (dmagList != null) { 
				int dmagListSize = dmagList.size();
				double [] dmagArray = new double[dmagListSize];
				for (int mm=0; mm< dmagListSize; mm++) {
					dmagArray[mm] = dmagList.get(mm);
				}
				//Calculate median...
				Arrays.sort(dmagArray);
				int result = dmagArray.length % 2;
				double median = 0;
				if (result == 0) {
					int rightNumber = dmagArray.length / 2;
					int leftNumber = rightNumber - 1;
					median = (dmagArray[rightNumber] + dmagArray[leftNumber]) / 2;
				} else {
					median = dmagArray[dmagArray.length / 2];
				}
				//Calculate mean and error in the mean...
				double sum = 0;
				double sum2 = 0.;
				double mean = -1000;
				double sdMean = -1000;
				for (int mm=0; mm< dmagListSize; mm++) {
					sum = sum + dmagArray[mm];
					sum2 = sum2 + dmagArray[mm]*dmagArray[mm];
				}
				int ntot = dmagArray.length;
				mean = sum / ntot;
				if (ntot > 1) {
					double var = sum2/ntot - mean*mean;
					double ratio = (double) ntot / (double) (ntot-1);
					double sampleVar = ratio*var;
					sdMean = Math.sqrt(sampleVar/ntot);
				}
				rro.setDmagListMean(mean);
				rro.setDmagListSDMean(sdMean);
				rro.setDmagListMedian(median);
				rro.setDmagListSize(ntot);
				if (verbose > 0 && kk % 100 == 0) {
					System.out.println(regionid1 + "\t" + regionid2 + "\t" + dmagList.size() + "\t" + mean + "\t"+  sdMean + "\t" + median);
				} else if (verbose > 1) {
					System.out.println(regionid1 + "\t" + regionid2 + "\t" + dmagList.size() + "\t" + mean + "\t"+  sdMean + "\t" + median);
				}
				
			}
		
		}	

		
		//Remove overlaps with fewer than 5 matched stars...
		if (verbose > 0) {
			System.out.println("regionRawOffsetsListSize before culling:  " + regionRawOffsetsListSize);
		}
		for (int kk=0; kk<regionRawOffsetsListSize; kk++) {
			int kkinv = regionRawOffsetsListSize-kk-1;
			RegionRawOffsets rro = (RegionRawOffsets) regionRawOffsetsList.get(kkinv);
			int regionid1 = rro.getRegionid1();
			int regionid2 = rro.getRegionid2();
			int dmagListSize = rro.getDmagListSize();
			double dmagListMean = rro.getDmagListMean();
			double dmagListSDMean = rro.getDmagListSDMean();
			double dmagListMedian = rro.getDmagListMedian();
			if (dmagListSize < 5) {
				regionRawOffsetsList.remove(kkinv);
			}
		}				
		regionRawOffsetsListSize = regionRawOffsetsList.size();
		if (verbose > 0) {
			System.out.println("regionRawOffsetsListSize after culling:  " + regionRawOffsetsListSize);
		}
		
		if (verbose > 0) {
			for (int ii=0; ii<regionidListSize; ii++) {
				int region = (Integer) regionidList.get(ii);
				System.out.println(ii + " " + region + " " + regionQualityMap.get((Integer) region));
			}
		}
		
		// Initialize a bunch of stuff that sits in the iteration loop
		double chi2 = -1;
		double rms  = -1;
		int dof     = -1;
		
		double[] array1d = new double[regionidListSize];
		double[][] array2d = new double[regionidListSize][regionidListSize];
		double[][] dmagArray2d = new double[regionidListSize][regionidListSize];
		int[][] overlapArray2d = new int[regionidListSize][regionidListSize];
		
		DoubleMatrix2D AA = null;
		DoubleMatrix2D BB = null;
		DoubleMatrix2D XX = null;
		DoubleMatrix2D AAinv = null;
		DoubleMatrix2D II = null;
		
		//Initialize arrays...
		for (int ii=0; ii<regionidListSize; ii++) {
			array1d[ii] = 0.;
			for (int jj=0; jj<regionidListSize; jj++) {
				array2d[ii][jj] = 0.;
				overlapArray2d[ii][jj] = 0;
				dmagArray2d[ii][jj] = 0.;
			}
		}
		
		//Populate arrays...
		for (int kk=0; kk<regionRawOffsetsListSize; kk++) {
			RegionRawOffsets rro = (RegionRawOffsets) regionRawOffsetsList.get(kk);
			int regionid1 = rro.getRegionid1();
			int regionid2 = rro.getRegionid2();
			double dmagListMedian = rro.getDmagListMedian();
			int index1 = regionidList.indexOf((Integer) regionid1);
			int index2 = regionidList.indexOf((Integer) regionid2);
			//We have reversed the signs of the AA matrix entries relative
			// to what is given in Glazebrook et al. (1994).  In the 
			// Glazebrook et al. convention, we would do the following:
			// array2d[index1][index1]--
			// array2d[index1][index1]--
			//array2d[index1][index2] = 1;
			//array2d[index2][index1] = 1;
			array2d[index1][index1]++;
			array2d[index2][index2]++;
			array2d[index1][index2] = -1;
			array2d[index2][index1] = -1;
			array1d[index1] = array1d[index1] + dmagListMedian;
			array1d[index2] = array1d[index2] - dmagListMedian;
			overlapArray2d[index1][index2] = 1;
			overlapArray2d[index2][index1] = 1;
			dmagArray2d[index1][index2] =  dmagListMedian;
			dmagArray2d[index2][index1] = -dmagListMedian;
		}		

		//Reset values for regions with calibrated zeropoints
		for (int ii=0; ii<regionidListSize; ii++) {
			int region = ((Integer) regionidList.get(ii)).intValue();
			int regionQuality = ((Integer) regionQualityMap.get((Integer) region)).intValue();
			if (regionQuality == 1) {
				for (int jj=0; jj<regionidListSize; jj++) {
					array2d[ii][jj] = 0.;
				}
				array2d[ii][ii] = 1;
				array1d[ii] = 0.;
			}
		}
				
		//Now convert the arrays into matrices.
		AA = new SparseDoubleMatrix2D(regionidListSize,regionidListSize);
		BB = new DenseDoubleMatrix2D(regionidListSize,1);
		for (int ii=0; ii<regionidListSize; ii++) {
			BB.set(ii,0,array1d[ii]);
			for (int jj=0; jj<regionidListSize; jj++) {
				AA.set(ii,jj,array2d[ii][jj]);
			}
		}
		
		//QA:  create and output a FITS image file of the matrix AA.
		if (verbose > 0) {
			if (regionidListSize < 100) {
				System.out.println("Matrix AA: \n" + AA.toString());
				System.out.println("");
			}
			Fits f = new Fits();
			double[][] dimg = new double[regionidListSize][regionidListSize];
			for (int i=0; i<regionidListSize; i++) {
				for (int j=0; j<regionidListSize; j++) {
					dimg[i][j] = AA.get(i,j);
				}
			}				
			f.addHDU(Fits.makeHDU(dimg));
			BufferedFile bf = new BufferedFile("matrixAA.fits", "rw");
			f.write(bf);
			bf.flush();
			bf.close();
		}

		//QA:  create and output a FITS image file of the matrix BB.
		if (verbose > 0) {
			if (regionidListSize < 100) {
				System.out.println("Matrix BB: \n" + BB.toString());
				System.out.println("");
			}
			Fits f = new Fits();
			double[][] dimg = new double[regionidListSize][1];
			for (int i=0; i<regionidListSize; i++) {
				dimg[i][0] = BB.get(i,0);
			}
			f.addHDU(Fits.makeHDU(dimg));
			BufferedFile bf = new BufferedFile("matrixBB.fits", "rw");
			f.write(bf);
			bf.flush();
			bf.close();
		}
		
		//Perform inversion of matrix AA.
		if (verbose > 0) {
			System.out.println("start " + regionidListSize + "x" + regionidListSize + " matrix inversion");
			System.out.println("");
		}
		Algebra alg = new Algebra();
		Date date1 = new Date();
		AAinv = alg.inverse(AA);
		Date date2 = new Date();
		double timeSec = (date2.getTime()-date1.getTime())/1000.;
		if (verbose > 0) {
			DecimalFormat formatter = new DecimalFormat("0.##");
			System.out.println(regionidListSize + "x" + regionidListSize + " matrix inversion finished in " + formatter.format(timeSec) + " sec\n");
		}

		//QA:  create and output FITS image files of the inverted matrix AAinv 
		//     and of the identity matrix II found by multiplying AA by AAinv.		
		if (verbose > 0) {
			II = alg.mult(AA,AAinv);
			if (regionidListSize < 100) {
				System.out.println("Matrix AAinv: \n" + AAinv.toString());
				System.out.println("");
				System.out.println("Identity Matrix: \n" + II.toString());
				System.out.println("");
			}
			Fits f = new Fits();
			double[][] dimg = new double[regionidListSize][regionidListSize];
			for (int i=0; i<regionidListSize; i++) {
				for (int j=0; j<regionidListSize; j++) {
					dimg[i][j] = AAinv.get(i,j);
				}
			}
			f.addHDU(Fits.makeHDU(dimg));	
			BufferedFile bf = new BufferedFile("matrixAAinv.fits", "rw");
			f.write(bf);
			bf.flush();
			bf.close();

			for (int i=0; i<regionidListSize; i++) {
				for (int j=0; j<regionidListSize; j++) {
					dimg[i][j] = II.get(i,j);
				}
			}
			f.addHDU(Fits.makeHDU(dimg));		
			bf = new BufferedFile("matrixII.fits", "rw");
			f.write(bf);
			bf.flush();
			bf.close();

		}

		// Multiply AAinv by BB to obtain the fit values for the zeropoints of each region.
		XX = alg.mult(AAinv,BB);
		
		//QA:  create and output a FITS image file of the matrix XX.
		if (verbose > 0) {
			System.out.println("Matrix XX: \n" + XX.toString());
			System.out.println("");
			Fits f = new Fits();
			double[][] dimg = new double[regionidListSize][1];
			for (int ii=0; ii<regionidListSize; ii++) {
				dimg[ii][0] = XX.get(ii,0);
			}
			f.addHDU(Fits.makeHDU(dimg));
			BufferedFile bf = new BufferedFile("matrixXX.fits", "rw");
			f.write(bf);
			bf.flush();
			bf.close();
		}
		
		//Estimate rms of solution and rms of each region's zeropoint.
		//Output each region's zeropoint and zeropoint rms to outputFileName.
		File outputFile = new File(outputFileName);
		FileWriter writer = new FileWriter(outputFile);
		writer.write("#To correct the mags in each regionid, *subtract* the following zeropoints.\n");
		writer.write("#regionid     zeropoint     rms \n");
		double rmsAll = -1.;
		double sumAll = 0.;
		double nAll   = 0.;
		for (int ii=0; ii<regionidListSize; ii++) {
			int regionid1 = Integer.parseInt(regionidList.get(ii).toString());
			double zp1 = XX.get(ii,0);
			double rms1 = -1.;
			double sum1 = 0.;
			int n = 0;
			for (int jj=0; jj<regionidListSize; jj++) {
				if (ii != jj) {
					if (overlapArray2d[ii][jj] == 1) {
						double zp2 = XX.get(jj,0);
						double dmagListMedian = dmagArray2d[ii][jj];
						//Under the Glazebrook et al. (1994) sign conventions, we would
						// subtract zp2 from zp1, not the reverse.
						sum1 = sum1 + (dmagListMedian-(zp1-zp2))*(dmagListMedian-(zp1-zp2));
						n++;
						sumAll = sumAll + (dmagListMedian-(zp1-zp2))*(dmagListMedian-(zp1-zp2));
						nAll++;
					}
				}
			}	
			if (n > 0) {
				rms1 = sum1/n;
				rms1 = Math.sqrt(rms1);
			}
			//This formatted input requireds java 1.5 or higher...
			String outputLine = String.format("%1$-10d   %2$8.3f   %3$8.3f \n", regionid1, zp1, rms1);
			if (verbose > 0) {
				System.out.print(outputLine);
			}
			writer.write(outputLine);
		}
		if (nAll > 0) {
			rmsAll = sumAll/nAll;
			rmsAll = Math.sqrt(rmsAll);
		}
		String outputLine = String.format("rms of solution: %1$8.3f \n", rmsAll);
		if (verbose > 0) {
			System.out.println(outputLine);
		}
		outputLine = "#\n#" + outputLine;
		writer.write(outputLine);
		writer.close();

		if (verbose > 0) {
			System.out.println("\nResults have been written to the file " + outputFileName + "\n");
		}

		if (verbose > 0) {
			System.out.println("End Time: \t" + (new Date()).toString());
			System.out.println("\nThat's all, folks!\n");
		}     

		return;
		
	}
	
	public String getStarMatchFileName() {
		return starMatchFileName;
	}
	
	public void setStarMatchFileName(String starMatchFileName) {
		this.starMatchFileName = starMatchFileName;
	}
	
	public String getOutputFileName() {
		return outputFileName;
	}
		
		public void setOutputFileName(String outputFileName) {
			this.outputFileName = outputFileName;
		}
		
		public int getVerbose() {
			return verbose;
		}
		
		public void setVerbose(int verbose) {
			this.verbose = verbose;
		}
		
		public Date getDate() {
			return date;
		}
		
		public void setDate(Date date) {
			this.date = date;
		}
		
		public String[] getFilterList() {
			return filterList;
		}
		
		public void setFilterList(String[] filterList) {
			this.filterList = filterList;
		}
		
	}
