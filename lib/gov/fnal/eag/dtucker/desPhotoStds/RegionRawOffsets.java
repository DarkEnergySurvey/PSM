package gov.fnal.eag.dtucker.desPhotoStds;

import cern.colt.list.DoubleArrayList;

public class RegionRawOffsets {
	
	private int regionid1 = -1;
	private int regionid2 = -1;
	private DoubleArrayList dmagList = null;
	private double dmagListMean = -1000.;
	private double dmagListSDMean = -1000.;
	private double dmagListMedian = -1000.;
	private int dmagListSize = -1;
	
	public int getRegionid1() {
		return regionid1;
	}
	
	public void setRegionid1(int imageid1) {
		this.regionid1 = imageid1;
	}
	
	public int getRegionid2() {
		return regionid2;
	}
	
	public void setRegionid2(int imageid2) {
		this.regionid2 = imageid2;
	}

	public DoubleArrayList getDmagList() {
		return dmagList;
	}

	public void setDmagList(DoubleArrayList dmagList) {
		this.dmagList = dmagList;
	}
	
	public double getDmagListMean() {
		return dmagListMean;
	}

	public void setDmagListMean(double dmagListMean) {
		this.dmagListMean = dmagListMean;
	}

	public double getDmagListMedian() {
		return dmagListMedian;
	}

	public void setDmagListMedian(double dmagListMedian) {
		this.dmagListMedian = dmagListMedian;
	}

	public int getDmagListSize() {
		return dmagListSize;
	}

	public void setDmagListSize(int dmagListSize) {
		this.dmagListSize = dmagListSize;
	}

	public double getDmagListSDMean() {
		return dmagListSDMean;
	}

	public void setDmagListSDMean(double dmagListSDMean) {
		this.dmagListSDMean = dmagListSDMean;
	}

}
