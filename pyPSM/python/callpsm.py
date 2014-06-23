#!/usr/bin/env python
"""
Main pyPSM script.


Examples:

callpsm --help

callpsm.py -s db-destest --inputCatListFile psm-20140202-g-r47p01.inputcat.csv --outputSolutionFile psm-20140202-g-r47p01.solution.fits --outputCatListFile psm-20140202-g-r47p01.provenance.csv --stdcat standard_stars_v6.csv --ksolve --bsolve --verbose 2

"""

import sys
import os
import time
import datetime

def main():

	import argparse
	import time

	"""Create command line arguments"""
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('--section','-s',default='db-destest',
			    help='section in the .desservices file w/ DB connection information')
	parser.add_argument('--inputCatListFile', help='CSV file containing input list of catalog file names (plus other relevant data)', default='psmInputCatList.csv')
	parser.add_argument('--outputSolutionFile', help='FITS table containing results of PSM fit', default='psmFit.fits')
	parser.add_argument('--outputCatListFile', help='CSV file containing list of catalog files actually used in PSM solution (for purposes of tracking provenance)', default='psmOutputCatList.fits')
	parser.add_argument('--stdcat', help='standard star catalog', default='standard_stars_v6.csv')
	parser.add_argument('--bsolve',help='solve for b terms', default=False, action='store_true')
	parser.add_argument('--ksolve',help='solve for k term', default=False, action='store_true')
	parser.add_argument('--verbose', help='verbosity level of output to screen', default=0)

        args = parser.parse_args()

        callpsm(args)


def callpsm(args):
   
        import os

	print args
	print args.nite
	print args.section
	print args.band

	queryTextFile = "obsquery.%s.%s.txt" % (args.nite, args.band)

	queryText = "\
SELECT e.expnum, s.filename, \n \
       s.object_number, s.x_image, s.y_image, s.radeg, s.decdeg, s.flux_psf, s.fluxerr_psf, \n \
       3600.*s.fwhm_world as fwhm_arcsec, s.class_star, s.spread_model, s.spreaderr_model, s.flags, \n \
       e.nite, e.object, e.band, i.ccdnum, i.airmass,  \n \
       e.mjd_obs, e.exptime, i.skybrite, i.skysigma, \n \
       i.elliptic as image_ellipt, 0.27*i.fwhm as image_fwhm_arcsec, \n \
       i.saturate as image_sat_level, c.filetype, i.crpix1, i.crpix2, i.naxis2, i.naxis2, rasicam.* \n \
FROM EXPOSURE e, IMAGE i, WDF, CATALOG c, SE_OBJECT s, gruendl.rasicam_decam@desoper rasicam \n \
WHERE e.expnum=i.expnum AND  \n \
       i.filename=wdf.parent_name AND \n \
       wdf.child_name=c.filename AND \n \
       c.filename=s.filename AND \n \
       e.filename=rasicam.exposurename AND \n \
       c.filetype='cat_finalcut' AND \n \
       s.reqnum=64 AND \n \
       ( (s.class_star > 0.8) OR ((s.spread_model + 3.*spreaderr_model) between -0.003 AND 0.003) ) AND \n \
       s.flux_psf > 2000. AND \n \
       s.flags < 3 AND \n \
       (UPPER(rasicam.source)='HEADER') AND (UPPER(rasicam.gskyphot)='T') AND \n \
       s.filename in (select distinct filename from se_object where nite=%s and reqnum = 64 and band='%s') \n \
ORDER BY s.radeg \n" % (args.nite,args.band)

	print queryText

	os.system("/bin/rm -f "+queryTextFile)
	ofd=open(queryTextFile,'w')
	ofd.write(queryText)
	ofd.close()


     ##cmd = "rm -f twomass-matchStarsRaDec2.cat"
     #cmd = "rm -f %s" % (twomasscat)
     #status = os.system(cmd)
     ##cmd = "%s -c %s %s -b 100.0 -smj -m %i >> twomass-matchStarsRaDec2.cat" % (find2mass,raFieldHMS,decFieldDMS,matchNStar)
     #cmd = "%s -c %s %s -b 100.0 -smj -m %i >> %s" % (find2mass,raFieldHMS,decFieldDMS,matchNStar,twomasscat)
     #print "Running %s" % cmd
     #status = os.system(cmd)
     #print status
     #if status != 0:
     #  print "find2mass failed on reference image %s" % image
     #  sys.exit(1)
     #else:
     #  print "find2mass ran successfully on reference image %s" % image

   # Output the results of fit...
#   os.system("/bin/rm -f "+inmatches+'.'+fitband+'.tmp*')
#   ofd=open(outak,'w')
#   ofd.write('niter_'+band+' '+str('%d'%int(niter))+'\n')
#   ofd.write('thresholdit_'+band+' '+str('%.4f'%float(thresholdit))+'\n')
#   ofd.write('n_'+band+' '+str('%d'%int(ninfit))+'\n')
#   ofd.write('dof_'+band+' '+str('%d'%int(dof))+'\n')
#   ofd.write('rms_'+band+' '+str('%.4f'%float(sigma))+'\n')
#   ofd.write('chisq_'+band+' '+str('%.4f'%float(chisq))+'\n')
#   ofd.write('k_'+band+' '+str('%.3f'%float(XX[0]))+' +/- '+str('%.3f'%float(errors[0]))+'\n')
#   for i in range(1,63):
#      outputLine = 'a_%s %2d %.3f +/- %.3f \n' % (band, i, XX[i]-25., errors[i])
#      ofd.write(outputLine)
#      if verbose > 1: print outputLine,
#   #endfor
#   for i in range(1,63):
#      outputLine = 'b_%s %2d %.3f +/- %.3f \n' % (band, i, XX[nccd+i], errors[nccd+i])
#      ofd.write(outputLine)
#      if verbose > 1: print outputLine,
#   #endfor
#   ofd.close()

#	query.py --section db-destest --format csv --header + < obsquery.$band.txt | awk 'NR>1' > obsquery.$band.out
#	echo "expnum,ccd,exptime,airmass,amp,ra,dec,mag,magerr,class_star,flags,band" > mystars.$band.csv
#	# Note:  the values of fluxerr_psf seem abnormally low; we will calculate the magerr_psf directly from the flux and a gain of 4.8e-/ADU...
#        #awk -F, 'BEGIN{OFS=","} NR>1 && $8>0. {print $1,$18,$21,$19,"-1",$6,$7,-2.5*0.43429448190325176*log($8)+25.,2.5*0.43429448190325176*$9/$8,$11,$14,$17}' obsquery.$band.out >> mystars.$band.csv
#	awk -F, 'BEGIN{OFS=","} NR>1 && $8>0. {print $1,$18,$21,$19,"-1",$6,$7,-2.5*0.43429448190325176*log($8)+25.,2.5*0.43429448190325176/sqrt(4.8*$8),$11,$14,$17}' obsquery.$band.out >> mystars.$band.csv
#	matchsort.py eqstds.csv mystars.$band.csv matchemup.$band.csv
#	psm.py matchemup.$band.csv solve.$band.txt $bandIndex 3 0.1 --bsolve --ksolve --verbose 3



#	echo "running QA..."
#	psmQA.py matchemup.$band.csv.$band.res.csv
#	mv PSM_QA_res_vs_airmass.png PSM_QA_res_vs_airmass.$band.png
#	mv PSM_QA_res_vs_ccd.png PSM_QA_res_vs_ccd.$band.png
#	mv PSM_QA_res_vs_color.png PSM_QA_res_vs_color.$band.png
#	mv PSM_QA_res_vs_expnum.png PSM_QA_res_vs_expnum.$band.png
#	mv PSM_QA_res_vs_mag.png PSM_QA_res_vs_mag.$band.png
	
#	echo "the solution for" $band "can be found in "matchemup.csv.$band.fit



if __name__ == "__main__":
	main()
