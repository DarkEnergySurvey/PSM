package gov.fnal.eag.dtucker.desPhotoStds;

public class MatchedStdStarYear1 {
	
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
	private double stdzY        = -1000.0;
	private double mjd          = -1000.0;
	private double stdRA        = -1000.0;
	private double stdDec       = -1000.0;
	private double fieldRA      = -1000.0;
	private double fieldDec     = -1000.0; 
	private long image_id       = -1000;
	private long exposure_id    = -1000;
	private double x_image      = -1000.0;
	private double y_image      = -1000.0;
	private long object_id      = -1;

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
	public double getStdzY() {
		return stdzY;
	}
	public void setStdzY(double stdzY) {
		this.stdzY = stdzY;
	}
	public double getStdRA() {
		return stdRA;
	}
	public void setStdRA(double stdRA) {
		this.stdRA = stdRA;
	}
	public double getStdDec() {
		return stdDec;
	}
	public void setStdDec(double stdDec) {
		this.stdDec = stdDec;
	}
	public double getFieldRA() {
		return fieldRA;
	}
	public void setFieldRA(double fieldRA) {
		this.fieldRA = fieldRA;
	}
	public double getFieldDec() {
		return fieldDec;
	}
	public void setFieldDec(double fieldDec) {
		this.fieldDec = fieldDec;
	}
	public long getImage_id() {
		return image_id;
	}
	public void setImage_id(long image_id) {
		this.image_id = image_id;
	}
	public double getX_image() {
		return x_image;
	}
	public void setX_image(double x_image) {
		this.x_image = x_image;
	}
	public double getY_image() {
		return y_image;
	}
	public void setY_image(double y_image) {
		this.y_image = y_image;
	}
	public long getObject_id() {
		return object_id;
	}
	public void setObject_id(long object_id) {
		this.object_id = object_id;
	}
	public long getExposure_id() {
		return exposure_id;
	}
	public void setExposure_id(long exposure_id) {
		this.exposure_id = exposure_id;
	}
	
}
