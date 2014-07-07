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
    parser.add_argument('--nite', help='nite to be queried', default=20131002)
    parser.add_argument('--reqnum', help='request number to be queried', default=3, type=int)
    parser.add_argument('--attempt', help='processing attempt to be queried (***CURRENTLY NOT USED***)', default=1, type=int)
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
    queryText = """
        SELECT e.expnum, concat('/full/path/',s.filename) as filename,
               e.nite, e.object, e.band, i.ccdnum, i.airmass,
               e.mjd_obs, e.exptime, i.skybrite, i.skysigma,
               i.elliptic as image_ellipt, 0.27*i.fwhm as image_fwhm_arcsec,
               i.saturate as image_sat_level, c.filetype, i.crpix1, i.crpix2, i.naxis2, i.naxis2,
               rasicam.gskyphot, rasicam.gskyvar, rasicam.lskyphot, rasicam.lskyvar, rasicam.source
        FROM EXPOSURE e, IMAGE i, WDF, CATALOG c, gruendl.rasicam_decam@desoper rasicam,
             (select distinct filename from se_object where nite=%s and reqnum=%d and band='%s') s
        WHERE e.expnum=i.expnum AND
              i.filename=wdf.parent_name AND
              wdf.child_name=c.filename AND
              c.filename=s.filename AND
              e.filename=rasicam.exposurename AND
              c.filetype='cat_finalcut'
        ORDER BY s.filename""" % (args.nite,args.reqnum,args.band)
    
        #AND
        #   (UPPER(rasicam.source)='HEADER') AND (UPPER(rasicam.gskyphot)='T')
    
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
