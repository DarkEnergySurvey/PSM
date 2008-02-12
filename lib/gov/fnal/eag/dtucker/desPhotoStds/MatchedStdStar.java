package gov.fnal.eag.dtucker.desPhotoStds;

public class MatchedStdStar {
	
	//Instance variables
	private String stdStarName  = "";
	private String fieldName    = "";
	private int counter			 = -1000;
	private int ccd_number      = -1000;
	private double airmass      = -1000.0;
	private double deltamag     = -1000.0;
	private double deltamagerr  = -1000.0;
	private double stdmag       = -1000.0;
	private double stdug        = -1000.0;
	private double stdgr        = -1000.0;
	private double stdri        = -1000.0;
	private double stdiz        = -1000.0;
	private double mjd          = -1000.0;

	public double getAirmass() {
		return airmass;
	}
	public void setAirmass(double airmass) {
		this.airmass = airmass;
	}
	public int getCcd_number() {
		return ccd_number;
	}
	public void setCcd_number(int ccd_number) {
		this.ccd_number = ccd_number;
	}
	public double getDeltamag() {
		return deltamag;
	}
	public void setDeltamag(double deltamag) {
		this.deltamag = deltamag;
	}
	public double getMjd() {
		return mjd;
	}
	public void setMjd(double mjd) {
		this.mjd = mjd;
	}
	public double getStdmag() {
		return stdmag;
	}
	public void setStdmag(double stdmag) {
		this.stdmag = stdmag;
	}
	public double getDeltamagerr() {
		return deltamagerr;
	}
	public void setDeltamagerr(double deltamagerr) {
		this.deltamagerr = deltamagerr;
	}
	public double getStdgr() {
		return stdgr;
	}
	public void setStdgr(double stdgr) {
		this.stdgr = stdgr;
	}
	public double getStdiz() {
		return stdiz;
	}
	public void setStdiz(double stdiz) {
		this.stdiz = stdiz;
	}
	public double getStdri() {
		return stdri;
	}
	public void setStdri(double stdri) {
		this.stdri = stdri;
	}
	public double getStdug() {
		return stdug;
	}
	public void setStdug(double stdug) {
		this.stdug = stdug;
	}
	public int getCounter() {
		return counter;
	}
	public void setCounter(int counter) {
		this.counter = counter;
	}
	public String getFieldName() {
		return fieldName;
	}
	public void setFieldName(String fieldName) {
		this.fieldName = fieldName;
	}
	public String getStdStarName() {
		return stdStarName;
	}
	public void setStdStarName(String stdStarName) {
		this.stdStarName = stdStarName;
	}
	
}
