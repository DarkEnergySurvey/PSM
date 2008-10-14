/**
 *  ColorTermCoeffs.java dtucker Sep 18, 2008 3:32:53 PM
 */
package gov.fnal.eag.dtucker.desPhotoStds;

import java.util.ArrayList;

/**
 * ColorTermCoeffs holds arrays describing the default values of the
 * color term ("b term") coefficients used in the photometric equations
 * @author dtucker
 * 
 */
public class ColorTermCoeffs {

	//Instance variables
	private ArrayList bccdidArrayList = new ArrayList();
	private ArrayList bdefaultArrayList = new ArrayList();
	private ArrayList bdefaultErrArrayList = new ArrayList();
	

	/**
	 * @return Returns the bccdidArrayList.
	 */
	public ArrayList getBccdidArrayList() {
		return bccdidArrayList;
	}

	/**
	 * @param bccdidArrayList The bccdidArrayList to set.
	 */
	public void setBccdidArrayList(ArrayList ccdArrayList) {
		this.bccdidArrayList = ccdArrayList;
	}	

	/**
	 * @return Returns the bdefaultArrayList.
	 */
	public ArrayList getBdefaultArrayList() {
		return bdefaultArrayList;
	}

	/**
	 * @param bdefaultArrayList The bdefaultArrayList to set.
	 */
	public void setBdefaultArrayList(ArrayList bdefaultArrayList) {
		this.bdefaultArrayList = bdefaultArrayList;
	}

	/**
	 * @return Returns the bdefaultErrArrayList.
	 */
	public ArrayList getBdefaultErrArrayList() {
		return bdefaultErrArrayList;
	}

	/**
	 * @param bdefaultErrArrayList The bdefaultErrArrayList to set.
	 */
	public void setBdefaultErrArrayList(ArrayList bdefaultErrArrayList) {
		this.bdefaultErrArrayList = bdefaultErrArrayList;
	}

}
