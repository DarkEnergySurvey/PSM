#!/usr/bin/env python
"""
Description:

Main pyPSM script.


Examples:

callpsm --help

callpsm.py -s db-destest --inputCatListFile psm-20140202-g-r47p01.inputcat.csv --outputSolutionFile psm-20140202-g-r47p01.solution.fits --outputCatListFile psm-20140202-g-r47p01.provenance.csv --stdcat standard_stars_v6.csv --ksolve --bsolve --verbose 2

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
    parser.add_argument('--band', help='band to be fit (TO BE DEPRECATED)', default='g', choices=['u','g','r','i','z','Y'])
    parser.add_argument('--magType', help='mag type to use (mag_psf, mag_auto, mag_aper_8, ...)', default='mag_psf')
    parser.add_argument('--sex_mag_zeropoint', help='default sextractor zeropoint to use to convert fluxes to sextractor mags (mag_sex = -2.5log10(flux) + sex_mag_zeropoint)', type=float, default=25.0)
    parser.add_argument('--bsolve',help='a toggle to solve for the b term coefficients', default=False, action='store_true')
    parser.add_argument('--ksolve',help='a toggle to solve for the k term coefficient', default=False, action='store_true')
    parser.add_argument('--niter', help='the number of iterations for the outlier rejection', type=int, default=0)
    parser.add_argument('--thresholdit', help='the threshold (in mag) for the outlier rejection', type=float, default=0)
    parser.add_argument('--magLo', help='std star mag lower limit  to use in fit', type=float, default=15.0)
    parser.add_argument('--magHi', help='std star mag upper limit to use in fit', type=float, default=18.0)
    parser.add_argument('--exptimeLo', help='exposure time lower limit (in seconds) to use in fit', type=float, default=2.0)
    parser.add_argument('--exptimeHi', help='exposure time upper limit (in seconds) to use in fit', type=float, default=100.0)
    parser.add_argument('--ccdExcludeList', help='a comma-separated list of CCDNUMs to exclude from the fit', default='')
    parser.add_argument('--expnumExcludeList', help='a comma-separated list of EXPNUMs to exclude from the fit', default='')
    parser.add_argument('--ignoreRasicam',help='include this flag to ignore RASICAM sky condition indicators', default=False, action='store_true')
    parser.add_argument('--keepIntermediateFiles',help='include this flag to keep (non-essential) intermediate files after running script', default=False, action='store_true')
    parser.add_argument('--project', help='the project id (deprecated?)', default='OPS')
    parser.add_argument('--psmversion', help='the version of pyPSM being used (deprecated?)', default='pyPSM_v0.1')
    parser.add_argument('--verbose', help='verbosity level of output to screen', default=0)

    args = parser.parse_args()

    if args.verbose > 0: print args

    callpsm(args)


#---------------------------------------------------------------------------

def callpsm(args):

    import string

    print args

#---------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#---------------------------------------------------------------------------
