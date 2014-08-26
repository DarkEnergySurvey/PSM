#!/usr/bin/env python
"""
    To be replaced by Michelle's query...
    
    create-cataloglist-from-db.py outputs a CSV list of the names of the single-epoch catalogs
    and additional exposure- and image-based parameters that are to be fed to the psm match
    and fit steps downstream of this script...
    
    
    Example:
    
    create-cataloglist-from-db.py --help
    
    create-cataloglist-from-db.py -s db-destest --nite 20131002 --reqnum 3 --attempt 1 --band g --outputCatListFile psmcats-20131002-g-r03p01.list --verbose 1
    
    """

import os
import time
import datetime


##################################

def main():
    import argparse
    import time
    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--section','-s',default='db-destest',
                        help='section in the .desservices file w/ DB connection information')
    parser.add_argument('--nite', help='nite to be queried', default='20140819')
    parser.add_argument('--tag', help='tag string number to be queried', default='Y2N_FIRSTCUT')
    parser.add_argument('--band', help='band to be queried', default='g', choices=['u','g','r','i','z','Y'])
    parser.add_argument('--outputCatListFile', help='ASCII list file containing list of unique names of all catalog files (plus exposure and image-based info) returned by query', default='psmcats.list')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
                        
    args = parser.parse_args()
                        
    if args.verbose > 0: print args
                        
    callcreatecataloglistfromdb(args)


##################################

def callcreatecataloglistfromdb(args):
    
    import coreutils.desdbi
    import cx_Oracle
    import os
    import csv
    import time
    
    #queryText = """SELECT DISTINCT filename FROM se_object WHERE nite=%s AND reqnum = 64 and BAND='%s' ORDER BY filename""" % (args.nite,args.band)
    
    # We add "/full/path/" to the catalog file name as a test for downstream
    # (Michelle gave a heads up that catalog names might include the path as well
    # as the base name).
    # NOTE:  RASICAM cuts (if any) are now performed downstream (in extract-catalog-from-db.py).
    #        The RASICAM info is grabbed via this query, but no RASICAM cuts now appear in the
    #        the WHERE clause of this query.
    queryText = """SELECT e.expnum, s.filename, e.nite, e.object, e.band, i.ccdnum, i.airmass,
 	e.mjd_obs, e.exptime, i.skybrite, i.skysigma, i.elliptic as image_ellipt, 
	0.27*i.fwhm as image_fwhm_arcsec, i.saturate as image_sat_level, c.filetype, 
	i.crpix1, i.crpix2, i.naxis1, i.naxis2 ,rasicam.gskyphot, rasicam.gskyvar, 
	rasicam.lskyphot, rasicam.lskyvar, rasicam.source 
	FROM PROD.EXPOSURE e left outer join PROD.RASICAM_DECAM rasicam on e.expnum = rasicam.expnum, 
	PROD.IMAGE i, PROD.WDF, PROD.CATALOG c, (select c.filename from PROD.catalog c,
 	PROD.ops_proctag optag, PROD.exposure e where e.expnum=c.expnum and 
	c.filename like concat(optag.unitname,'%s') and c.filetype='cat_finalcut' and 
	optag.tag='%s' and c.nite='%s' and c.band='%s' and 
	e.program='photom-std-field' order by c.expnum, c.ccdnum) s WHERE e.expnum=i.expnum AND 
	i.filename=wdf.parent_name AND wdf.child_name=c.filename AND c.filename=s.filename ORDER BY 
        s.filename""" % ('%',args.tag,args.nite,args.band)

    if args.verbose > 0: print queryText
    
    print 'Running query...'
    dbh = coreutils.DesDbi(os.path.join(os.getenv("HOME"),".desservices.ini"),args.section)
    t0=time.time()
    cur = dbh.cursor()
    cur.execute(queryText)
    
    if args.verbose > 0: print "query took", time.time()-t0, "seconds"
    
    print "Outputting query results to file", args.outputCatListFile
    with open(args.outputCatListFile,'w') as csvFile:
        writer = csv.writer(csvFile,delimiter=',',  quotechar='|',
                            lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
        #First output header...
        hdr = [col[0] for col in cur.description]
        writer.writerow(hdr)
                            
        #Then output contents...
        for line in cur:
            writer.writerow(line)
    
    dbh.close()


##################################


if __name__ == "__main__":
    main()


##################################
