#!/usr/bin/env python

# Author:    Douglas Tucker
# Date:      11 Jul 2014

"""
Description:

Main pyPSM script.

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* WARNING:  callpsm.py assumes that the name of the inputCatFileList has the form "psmcats-%8d-%1s-r%2dp%2d.list" *
*           (e.g., "psmcats-20131002-g-r03p01.list") in order to parse names for intermediate files correctly     *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

Examples:

callpsm --help

callpsm.py -s db-destest --inputCatListFile psmcats-20131002-g-r03p01.list --outputSolutionFile psmResults-20131002-g-r03p01.fits --outputCatListFile psmcats-20131002-g-r03p01.provenance.list --stdcat standard_stars_all_id6_pypsmFormat.csv --niter 3 --thresholdit 0.1 --ksolve --bsolve --ccdExcludeList 61,2,31 --expnumExcludeList 240688,240724,240774 --verbose 2

"""

import sys
import os
import time
import datetime

#---------------------------------------------------------------------------

def main():

    import argparse

    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--section','-s',default='db-destest',help='section in the .desservices file w/ DB connection information')
    parser.add_argument('--inputCatListFile', help='CSV file containing input list of catalog file names (plus other relevant data)', default='psmInputCatList.csv')
    parser.add_argument('--outputSolutionFile', help='FITS table containing results of PSM fit', default='psmFit.fits')
    parser.add_argument('--outputCatListFile', help='CSV file containing list of catalog files actually used in PSM solution (for purposes of tracking provenance)', default='psmOutputCatList.fits')
    parser.add_argument('--stdcat', help='standard star catalog', default='standard_stars_pypsmFormat.csv')
    parser.add_argument('--magType', help='mag type to use (mag_psf, mag_auto, mag_aper_8, ...)', default='mag_psf')
    parser.add_argument('--sex_mag_zeropoint', help='default sextractor zeropoint to use to convert fluxes to sextractor mags (mag_sex = -2.5log10(flux) + sex_mag_zeropoint)', type=float, default=25.0)
    parser.add_argument('--bsolve',help='a toggle to solve for the b term coefficients', default=False, action='store_true')
    parser.add_argument('--ksolve',help='a toggle to solve for the k term coefficient', default=False, action='store_true')
    parser.add_argument('--niter', help='the number of iterations for the outlier rejection', type=int, default=0)
    parser.add_argument('--thresholdit', help='the threshold (in mag) for the outlier rejection', type=float, default=0)
    parser.add_argument('--magLo', help='std star mag lower limit  to use in fit', type=float, default=15.0)
    parser.add_argument('--magHi', help='std star mag upper limit  to use in fit', type=float, default=18.0)
    parser.add_argument('--exptimeLo', help='exposure time lower limit (in seconds) to use in fit', type=float, default=2.0)
    parser.add_argument('--exptimeHi', help='exposure time upper limit (in seconds) to use in fit', type=float, default=100.0)
    parser.add_argument('--ccdExcludeList', help='a comma-separated list of CCDNUMs to exclude from the fit', default='')
    parser.add_argument('--expnumExcludeList', help='a comma-separated list of EXPNUMs to exclude from the fit', default='')
    parser.add_argument('--ignoreRasicam',help='include this flag to ignore RASICAM sky condition indicators', default=False, action='store_true')
    parser.add_argument('--keepIntermediateFiles',help='include this flag to keep (non-essential) intermediate files after running script', default=False, action='store_true')
    parser.add_argument('--project', help='the project id (deprecated?)', default='OPS')
    parser.add_argument('--psmversion', help='the version of pyPSM being used (deprecated?)', default='pyPSM_v0.1')
    parser.add_argument('--verbose', help='verbosity level of output to screen', type=int, default=0)

    args = parser.parse_args()

    if args.verbose > 0: print args

    callpsm(args)


#---------------------------------------------------------------------------

def callpsm(args):

    import numpy as np
    import string

    # parameters from original callpsm.py args list...
    section = args.section
    inputCatListFile = args.inputCatListFile
    outputSolutionFile  = args.outputSolutionFile
    outputCatListFile = args.outputCatListFile
    stdcat = args.stdcat
    magType = args.magType
    sex_mag_zeropoint = args.sex_mag_zeropoint
    bsolve = args.bsolve
    ksolve = args.ksolve
    niter = args.niter
    thresholdit = args.thresholdit
    magLo = args.magLo
    magHi = args.magHi
    exptimeLo = args.exptimeLo
    exptimeHi = args.exptimeHi
    ccdExcludeList = args.ccdExcludeList
    expnumExcludeList = args.expnumExcludeList
    ignoreRasicam = args.ignoreRasicam
    keepIntermediateFiles = args.keepIntermediateFiles
    project = args.project
    psmversion = args.psmversion
    verbose = args.verbose

    # Some quantities derived from the parameters from original callpsm.py args list...

    # There's got to be a better way to extract the year/band/reqnum/processid 
    # section of the inputCatListFile (e.g., something like "scanf"), but this
    # is a quick-and-dirty solution for now...
    baseName = inputCatListFile[8:25]

    # Likewise for the filter band name, which is currently needed by psm.py
    # (although psm.py should just extract it from the contents of the 
    # inputCatListFile table itself...).  Again, a quick-and-dirty solution
    # for now...
    band = inputCatListFile[17:18]

    # The following will be used below to determine how many -- if any -- expums/ccds are in the 
    #  expnumExcludeList/ccdExcludeList...
    expnumExcludeArray = np.fromstring(expnumExcludeList,dtype=int,sep=',')
    ccdExcludeArray = np.fromstring(ccdExcludeList,dtype=int,sep=',')

    print '\n\n\n'

    # First command:  extract-catalog-from-db.py
    outputObsFile = 'obsquery-%s.csv' % baseName
    cmd = """extract-catalog-from-db.py -s %s --inputCatListFile %s --outputObsFile %s --magType %s --sex_mag_zeropoint %f --verbose %d""" % (section, inputCatListFile, outputObsFile, magType, sex_mag_zeropoint, verbose)
    if ignoreRasicam:
        cmd = cmd+" --ignoreRasciam"
    if verbose > 0:
        print 'Running:  '+cmd
    status = os.system(cmd)
    if status != 0:
        print '* * * extract-catalog-from-db.py failed... exiting callpsm.py now! * * *'
        sys.exit(1)
    else:
        if verbose > 0: print '\n\n\n'


    # Second command:  matchsort.py
    inputObsCatFile = outputObsFile
    outputMatchFile = 'matched-%s.csv' % baseName
    cmd = """matchsort.py --inputStdStarCatFile %s --inputObsCatFile %s --outputMatchFile %s --verbose %d""" % (stdcat, inputObsCatFile, outputMatchFile, verbose)
    if verbose > 0:
        print 'Running:  '+cmd
    status = os.system(cmd)
    if status != 0:
        print '* * * matchsort.py failed... exiting callpsm.py now! * * *'
        sys.exit(1)
    else:
        if verbose > 0: print '\n\n\n'


    # Third command:  psm.py
    inputMatchFile = outputMatchFile
    outputResultsFITSFile = 'psmResults-%s.fits' % baseName
    outputResultsLogFile = 'psmResults-%s.log' % baseName
    outputResidualsFile = 'psmResiduals-%s.csv' % baseName
    outputCatsUsedFile = outputCatListFile
    cmd = """psm.py --inputMatchFile %s --outputResultsFITSFile %s --outputResultsLogFile %s --outputResidualsFile %s --outputCatsUsedFile %s --band %s --niter %d --thresholdit %f --exptimeLo %f --exptimeHi %f --magLo %f --magHi %f --project %s --psmversion %s --verbose %d""" % (inputMatchFile, outputResultsFITSFile, outputResultsLogFile, outputResidualsFile, outputCatsUsedFile, band, niter, thresholdit, exptimeLo, exptimeHi, magLo, magHi, project, psmversion, verbose)
    if bsolve:
        cmd = cmd+" --bsolve"
    if ksolve:
        cmd = cmd+" --ksolve"
    if expnumExcludeArray.size > 0:
        cmd = cmd+' --expnumExcludeList '+expnumExcludeList
    if ccdExcludeArray.size > 0:
        cmd = cmd+' --ccdExcludeList '+ccdExcludeList
    if verbose > 0:
        print 'Running:  '+cmd
    status = os.system(cmd)
    if status != 0:
        print '* * * psm.py failed... exiting callpsm.py now! * * *'
        sys.exit(1)
    else:
        if verbose > 0: print '\n\n\n'



    # Fourth command:  psmQA.py
    inputResidualsFile = outputResidualsFile
    outputFileBaseName = 'psmQA-%s' % baseName
    cmd = """psmQA.py --inputResidualsFile %s --outputFileBaseName %s --verbose %d""" % (inputResidualsFile, outputFileBaseName, verbose)
    if verbose > 0:
        print 'Running:  '+cmd
    status = os.system(cmd)
    if status != 0:
        print '* * * psmQA.py failed... exiting callpsm.py now! * * *'
        sys.exit(1)
    else:
        if verbose > 0: print '\n\n\n'


    # Fifth command:  cleanup intermediate files
    if keepIntermediateFiles:
        if verbose > 0:
            print 'Final Cleanup:  keeping intermediate files...'
            print '\n\n\n'
    else:
        if verbose > 0:
            print 'Final Cleanup:  deleting following intermediate files:  %s, $s' % (outputObsFile, outputMatchFile)
        cmd = """rm -f %s %s""" % (outputObsFile, outputMatchFile)
        status = os.system(cmd)
        if status != 0:
            print '* * * Final cleanup failed... exiting callpsm.py now! * * *'
            sys.exit(1)
        else:
            if verbose > 0: print '\n\n\n'
    #endif keepIntermediateFiles 


    # Final closeout...
    if verbose > 0:
        print 'callpsm.py finished successfully'

    sys.exit(0)


#---------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#---------------------------------------------------------------------------
