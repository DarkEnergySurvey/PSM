package gov.fnal.eag.dtucker.desPhotoStds;

public class MatchedStar {

	//Instance variables
	private int regionid1        = -1;
	private double regionRaCenDeg1    = -1000.;
	private double regionDecCenDeg1   = -1000.;
	private int regionQuality1   = -1;
	private int starid1         = -1;
	private double raDeg1       = -1000.;
	private double decDeg1      = -1000.;
	private double mag1         = -1000.;
	private double magErr1      = -1000.;
	
	private int regionid2        = -1;
	private double regionRaCenDeg2    = -1000.;
	private double regionDecCenDeg2   = -1000.;
	private int regionQuality2   = -1;
	private int starid2         = -1;
	private double raDeg2       = -1000.;
	private double decDeg2      = -1000.;
	private double mag2         = -1000.;
	private double magErr2      = -1000.;
	
	private double sepArcsec12  = -1000.;

	public double getDecDeg1() {
		return decDeg1;
	}

	public void setDecDeg1(double decDeg1) {
		this.decDeg1 = decDeg1;
	}

	public double getDecDeg2() {
		return decDeg2;
	}

	public void setDecDeg2(double decDeg2) {
		this.decDeg2 = decDeg2;
	}

	public double getMag1() {
		return mag1;
	}

	public void setMag1(double mag1) {
		this.mag1 = mag1;
	}

	public double getMag2() {
		return mag2;
	}

	public void setMag2(double mag2) {
		this.mag2 = mag2;
	}

	public double getMagErr1() {
		return magErr1;
	}

	public void setMagErr1(double magErr1) {
		this.magErr1 = magErr1;
	}

	public double getMagErr2() {
		return magErr2;
	}

	public void setMagErr2(double magErr2) {
		this.magErr2 = magErr2;
	}

	public double getRaDeg1() {
		return raDeg1;
	}

	public void setRaDeg1(double raDeg1) {
		this.raDeg1 = raDeg1;
	}

	public double getRaDeg2() {
		return raDeg2;
	}

	public void setRaDeg2(double raDeg2) {
		this.raDeg2 = raDeg2;
	}

	public double getRegionDecCenDeg1() {
		return regionDecCenDeg1;
	}

	public void setRegionDecCenDeg1(double regionDecCenDeg1) {
		this.regionDecCenDeg1 = regionDecCenDeg1;
	}

	public double getRegionDecCenDeg2() {
		return regionDecCenDeg2;
	}

	public void setRegionDecCenDeg2(double regionDecCenDeg2) {
		this.regionDecCenDeg2 = regionDecCenDeg2;
	}

	public int getRegionid1() {
		return regionid1;
	}

	public void setRegionid1(int regionid1) {
		this.regionid1 = regionid1;
	}

	public int getRegionid2() {
		return regionid2;
	}

	public void setRegionid2(int regionid2) {
		this.regionid2 = regionid2;
	}

	public int getRegionQuality1() {
		return regionQuality1;
	}

	public void setRegionQuality1(int regionQuality1) {
		this.regionQuality1 = regionQuality1;
	}

	public int getRegionQuality2() {
		return regionQuality2;
	}

	public void setRegionQuality2(int regionQuality2) {
		this.regionQuality2 = regionQuality2;
	}

	public double getRegionRaCenDeg1() {
		return regionRaCenDeg1;
	}

	public void setRegionRaCenDeg1(double regionRaCenDeg1) {
		this.regionRaCenDeg1 = regionRaCenDeg1;
	}

	public double getRegionRaCenDeg2() {
		return regionRaCenDeg2;
	}

	public void setRegionRaCenDeg2(double regionRaCenDeg2) {
		this.regionRaCenDeg2 = regionRaCenDeg2;
	}

	public double getSepArcsec12() {
		return sepArcsec12;
	}

	public void setSepArcsec12(double sepArcsec12) {
		this.sepArcsec12 = sepArcsec12;
	}

	public int getStarid1() {
		return starid1;
	}

	public void setStarid1(int starid1) {
		this.starid1 = starid1;
	}

	public int getStarid2() {
		return starid2;
	}

	public void setStarid2(int starid2) {
		this.starid2 = starid2;
	}

	
}
