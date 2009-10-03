package gov.fnal.eag.dtucker.desPhotoStds;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;

public class FocalPlaneCoords {
	
	//Instance variables
	private int nAxis1 = 0;
	private int nAxis2 = 0;
	private double xCenter = 0.0;
	private double yCenter = 0.0;
		
	public int getNAxis1() {
		return nAxis1;
	}
	public void setNAxis1(int axis1) {
		nAxis1 = axis1;
	}
	public int getNAxis2() {
		return nAxis2;
	}
	public void setNAxis2(int axis2) {
		nAxis2 = axis2;
	}
	public double getXCenter() {
		return xCenter;
	}
	public void setXCenter(double center) {
		xCenter = center;
	}
	public double getYCenter() {
		return yCenter;
	}
	public void setYCenter(double center) {
		yCenter = center;
	}
	

	
}
