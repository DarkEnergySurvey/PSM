#!/usr/bin/env python
"""
extract-catalog-from-db.py


Examples:

extract-catalog-from-db.py --help

extract-catalog-from-db.py -s db-destest --inputCatListFile psmcats-20131002-g-r03p01.list --outputObsFile obsquery-20131002-g-r03p01.csv --verbose 1
    
"""

import sys
import os
import time

##################################

def main():
    import argparse
    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--section','-s',default='db-destest',
                        help='section in the .desservices file w/ DB connection information')
    parser.add_argument('--inputCatListFile', help='CSV file containing list of unique names of all catalog files (and exposure- and image-based info) returned by query', default='psmcats.list')
    parser.add_argument('--outputObsFile', help='CSV file containing the observed star data to be matched and fit', default='obsdata.csv')
    parser.add_argument('--ignoreRasicam',help='include this flag to ignore RASICAM sky condition indicators', default=False, action='store_true')
    parser.add_argument('--keepIntermediateFiles',help='include this flag to keep (non-essential) intermediate files after running script', default=False, action='store_true')
    parser.add_argument('--magType', help='mag type to use (mag_psf, mag_auto, mag_aper_8, ...)', default='mag_psf')
    parser.add_argument('--sex_mag_zeropoint', help='default sextractor zeropoint to use to convert fluxes to sextractor mags (mag_sex = -2.5log10(flux) + sex_mag_zeropoint)', type=float, default=25.0)
    parser.add_argument('--verbose', help='verbosity level of output to screen (0, 1, 2, ...)', type=int, default=0)
                        
    args = parser.parse_args()

    callextractcatalogsfromdb(args)


##################################

def callextractcatalogsfromdb(args):
    
    import coreutils.desdbi
    import cx_Oracle
    import csv
    import numpy as np
    import math
    import string

    if args.verbose > 0: print args
    
    print 'Reading in catalog list input file %s' % args.inputCatListFile
    
    # Check first if file exists...
    if not os.path.isfile(args.inputCatListFile):
        print '%s does not seem to exist... exiting now...' % args.inputCatListFile
        sys.exit(1)
    
    # Read in the file...
    data=np.genfromtxt(args.inputCatListFile,dtype=None,delimiter=',',names=True)
    
    # Identify band...
    bandList = data['BAND'].tolist()
    sortedUniqueBandList = sorted(set(bandList))
    if len(sortedUniqueBandList) > 1:
        print 'Multiple bands in %s:' % args.inputCatListFile
        for sortedUniqueBand in sortedUniqueBandList:
            print sortedUniqueBand,
        print 'pyPSM cannot currently handle multiple bands... exiting now...'
        sys.exit(1)
    else:
        band = bandList[0]
        if args.verbose > 1:  print 'band = %s' % band
    
    # Identify nite...
    niteList = data['NITE'].tolist()
    sortedUniqueNiteList = sorted(set(niteList))
    if len(sortedUniqueNiteList) > 1:
        print 'Multiple nites in %s:' % args.inputCatListFile
        for sortedUniqueNite in sortedUniqueNiteList:
            print sortedUniqueNite,
        print 'pyPSM cannot currently handle multiple nites... exiting now...'
        sys.exit(1)
    else:
        nite = niteList[0]
        if args.verbose > 1:  print 'nite = %s' % nite
    
    # Create an array index mask based on RASICAM sky conditions...
    if args.ignoreRasicam:
        # Create a mask the with the same number of records as data['GSKYPHOT'],
        # but with all values set to Boolean "True"...
        rasicamMask = np.array((data['GSKYPHOT'].size)*[True],bool)
    else:
        # Use RASICAM info to include only those catalogs taken under apparently photometric conditions...
        rasicamMask = (data['SOURCE']=='HEADER') & (data['GSKYPHOT']=='T')
    
    # Create list of catalog file basenames...
    #  First, strip off path names from numpy array data['FILENAME']...
    for i in range(data['FILENAME'].size):
        data['FILENAME'][i] = os.path.basename(data['FILENAME'][i])
    # Then, convert the data['FILENAME'] numpy array into a python list...
    #  but only for those entries for which the rasicamMask is set to "True"...
    catFilenameList = data['FILENAME'][np.where(rasicamMask)].tolist()
    
    
    # Create python dictionaries from columns in the catalog input file using
    # entries from catFilenameList as the key...
    
    expnumList        = data['EXPNUM'][np.where(rasicamMask)].tolist()
    objectList        = data['OBJECT'][np.where(rasicamMask)].tolist()
    ccdnumList        = data['CCDNUM'][np.where(rasicamMask)].tolist()
    airmassList       = data['AIRMASS'][np.where(rasicamMask)].tolist()
    mjdobsList        = data['MJD_OBS'][np.where(rasicamMask)].tolist()
    exptimeList       = data['EXPTIME'][np.where(rasicamMask)].tolist()
    skybriteList      = data['SKYBRITE'][np.where(rasicamMask)].tolist()
    skysigmaList      = data['SKYSIGMA'][np.where(rasicamMask)].tolist()
    imageelliptList   = data['IMAGE_ELLIPT'][np.where(rasicamMask)].tolist()
    imageFWHMList     = data['IMAGE_FWHM_ARCSEC'][np.where(rasicamMask)].tolist()
    imagesatlevelList = data['IMAGE_SAT_LEVEL'][np.where(rasicamMask)].tolist()
    filetypeList      = data['FILETYPE'][np.where(rasicamMask)].tolist()
    crpix1List        = data['CRPIX1'][np.where(rasicamMask)].tolist()
    crpix2List        = data['CRPIX2'][np.where(rasicamMask)].tolist()
    #Fix following line after further tests...  (done!)
    #naxis1List        = data['NAXIS2'][np.where(rasicamMask)].tolist()
    naxis1List        = data['NAXIS1'][np.where(rasicamMask)].tolist()
    naxis2List        = data['NAXIS2'][np.where(rasicamMask)].tolist()
    gskyphotList      = data['GSKYPHOT'][np.where(rasicamMask)].tolist()
    gskyvarList       = data['GSKYVAR'][np.where(rasicamMask)].tolist()
    lskyphotList      = data['LSKYPHOT'][np.where(rasicamMask)].tolist()
    lskyvarList       = data['LSKYVAR'][np.where(rasicamMask)].tolist()
    sourceList        = data['SOURCE'][np.where(rasicamMask)].tolist()
    
    expnumDict        = dict(zip(catFilenameList, expnumList))
    objectDict        = dict(zip(catFilenameList, objectList))
    ccdnumDict        = dict(zip(catFilenameList, ccdnumList))
    airmassDict       = dict(zip(catFilenameList, airmassList))
    mjdobsDict        = dict(zip(catFilenameList, mjdobsList))
    exptimeDict       = dict(zip(catFilenameList, exptimeList))
    skybriteDict      = dict(zip(catFilenameList, skybriteList))
    skysigmaDict      = dict(zip(catFilenameList, skysigmaList))
    imageelliptDict   = dict(zip(catFilenameList, imageelliptList))
    imageFWHMDict     = dict(zip(catFilenameList, imageFWHMList))
    imagesatlevelDict = dict(zip(catFilenameList, imagesatlevelList))
    filetypeDict      = dict(zip(catFilenameList, filetypeList))
    crpix1Dict        = dict(zip(catFilenameList, crpix1List))
    crpix2Dict        = dict(zip(catFilenameList, crpix2List))
    naxis1Dict        = dict(zip(catFilenameList, naxis1List))
    naxis2Dict        = dict(zip(catFilenameList, naxis2List))
    gskyphotDict      = dict(zip(catFilenameList, gskyphotList))
    gskyvarDict       = dict(zip(catFilenameList, gskyvarList))
    lskyphotDict      = dict(zip(catFilenameList, lskyphotList))
    lskyvarDict       = dict(zip(catFilenameList, lskyvarList))
    sourceDict        = dict(zip(catFilenameList, sourceList))
    
    
    
    # We want to load the catFileNameList into the OPM_FILENAME_GTT table...
    
    print 'Opening connection to database...'
    dbh = coreutils.DesDbi(os.path.join(os.getenv("HOME"),".desservices.ini"),args.section)
    
    print 'Loading %d catalog filenames to OPM_FILENAME_GTT table...' % len(catFilenameList)
    dbh.load_filename_gtt(catFilenameList)
    

    # Create query...
    
    magType = args.magType
    magType = magType.strip()
    fluxType = magType.replace('mag','flux')
    fluxerrType = magType.replace('mag','fluxerr')
    
    queryText = """
        SELECT s.filename,
               s.object_number, s.x_image, s.y_image, s.radeg, s.decdeg, s.%s, s.%s,
               3600.*s.fwhm_world as fwhm_arcsec, s.class_star, s.spread_model, s.spreaderr_model, s.flags
        FROM SE_OBJECT s, opm_filename_gtt g
        WHERE  s.filename=g.filename AND
               ( (s.class_star > 0.8) OR ((s.spread_model + 3.*spreaderr_model) between -0.003 AND 0.003) ) AND
               s.flux_psf > 2000. AND
               s.flags < 3
        ORDER BY s.radeg""" % (fluxType, fluxerrType)
    
    
    if args.verbose > 0: print queryText
    
    print 'Running query...'
    
    t0=time.time()
    cur = dbh.cursor()
    cur.execute(queryText)
    
    if args.verbose > 0: print "query took", time.time()-t0, "seconds"
    

    print 'Outputting query results to file', args.outputObsFile
    
    with open(args.outputObsFile,'w') as csvFile:
        writer = csv.writer(csvFile,delimiter=',',  quotechar='|',
                            lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
                            
        #First output (modified) header...
        hdr = [col[0] for col in cur.description]
        hdr.insert(0,'EXPNUM')
        hdr.insert(9,'MAG')
        hdr.insert(10,'MAGERR')
        hdr.insert(11,'ZEROPOINT')
        hdr.insert(12,'MAGTYPE')
        hdr.extend(['NITE','OBJECT','BAND','CCDNUM','AIRMASS','MJD_OBS','EXPTIME','SKYBRITE','SKYSIGMA','IMAGE_ELLIPT','IMAGE_FWHM_ARCSEC','IMAGE_SAT_LEVEL','FILETYPE','CRPIX1','CRPIX2','NAXIS1','NAXIS2','GSKYPHOT','GSKYVAR','LSKYPHOT','LSKYVAR','SOURCE'])
        if args.verbose > 3:  print hdr
        writer.writerow(hdr)
                            
        #Then output (modified) contents...
        for line in cur:

            # First change line from a tuple to a list,
            # so we can add to it...
            line = list(line)

            # grab the catalog filename and place it in the "zeroth" column...
            catFilename = line[0].strip()
            line.insert(0,expnumDict[catFilename])

            # calculate the Sextractor-like instrumental mag and magerr from the flux and fluxerr
            #  and insert them (and the mag_zeropoint) just after the flux and fluxerr columns...
            # after inserting the catFileName to position 0 of line in the previous step, 
            #  line[7] should contain the flux and line[8] the fluxerr...
            if line[7] > 0.:
                mag = -2.5*math.log10(line[7]) + args.sex_mag_zeropoint
                magerr = (2.5/math.log(10.))*(line[8]/line[7])
            else:
                mag = -9999.
                magerr = -9999.

            line.insert(9,mag)
            line.insert(10,magerr)
            line.insert(11,args.sex_mag_zeropoint)
            line.insert(12,magType)

            line.extend([nite,
                         objectDict[catFilename],
                         band,
                         ccdnumDict[catFilename],
                         airmassDict[catFilename],
                         mjdobsDict[catFilename],
                         exptimeDict[catFilename],
                         skybriteDict[catFilename],
                         skysigmaDict[catFilename],
                         imageelliptDict[catFilename],
                         imageFWHMDict[catFilename],
                         imagesatlevelDict[catFilename],
                         filetypeDict[catFilename],
                         crpix1Dict[catFilename],
                         crpix2Dict[catFilename],
                         naxis1Dict[catFilename],
                         naxis2Dict[catFilename],
                         gskyphotDict[catFilename],
                         gskyvarDict[catFilename],
                         lskyphotDict[catFilename],
                         lskyvarDict[catFilename],
                         sourceDict[catFilename]
                         ])
            if args.verbose > 3:  print line
            writer.writerow(line)
    
    dbh.close()



##################################


if __name__ == "__main__":
    main()


##################################

