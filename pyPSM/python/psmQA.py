#!/usr/bin/env python

# Author:    Douglas Tucker
# Date:      21 May 2013
# Updated:   21 May 2013
# Updated:   11 Jul 2014

"""
Description:

This file defines methods for creating QA plots from the residuals
residuals CSV file output by psm.

Examples:

psmQA.py --help

psmQA.py --inputResidualsFile psmResiduals-20131002-g-r03p01.csv --outputFileBaseName psmQA-20131002-g-r03p01

"""

import numpy
import sys
import math
import os
import getopt
import csv
import matplotlib.pyplot as plt


#---------------------------------------------------------------------------

def main():
    import argparse
    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputResidualsFile', help='CSV file containing the residuals from psm.py', default='psmResiduals.csv')
    parser.add_argument('--outputFileBaseName', help='Base name for the QA plots that will be output)', default='psmQA')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0, 1, 2, ...)', type=int, default=0)

    args = parser.parse_args()

    if args.verbose > 0: print args

    psmQA(args)


#---------------------------------------------------------------------------

def psmQA(args):
    
    psmQATableFile = args.inputResidualsFile
    outBaseName = args.outputFileBaseName
    verbose = args.verbose

    #Read in csv file...
    csv_file_object = csv.reader(open(psmQATableFile, 'rb'))
    # Use this to skip cvs header...
    header = csv_file_object.next()
    data=[]
    for row in csv_file_object:
        data.append(row)
    #endfor
    
    # Convert from a list to an array.
    # (Be aware that each item is currently a string in this format)...
    data = numpy.array(data)
    
    # Create 1d arrays from the different columns of the CSV file,
    # converting to the appropriate format as necessary...
    res      = data[0::,0].astype(numpy.double)
    airmass  = data[0::,1].astype(numpy.double)
    magstd   = data[0::,2].astype(numpy.double)
    colorstd = data[0::,3].astype(numpy.double)
    ccd      = data[0::,4].astype(numpy.int)
    expnum   = data[0::,5].astype(numpy.int)
    mjdobs   = data[0::,6].astype(numpy.double)
    
    # Output QA plots...

    outputFile = '%s.res_vs_airmass.png' % (outBaseName)
    if verbose > 1:  print 'Outputting '+outputFile
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('residuals vs airmass')
    ax.set_xlabel('airmass',fontsize=12)
    ax.set_ylabel('residuals [mag]',fontsize=12)
    ax.plot(airmass, res, 'b.')
    ax.grid(True)
    plt.savefig(outputFile, format='png')
    
    outputFile = '%s.res_vs_mag.png' % (outBaseName)
    if verbose > 1:  print 'Outputting '+outputFile
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('residuals vs std mag')
    ax.set_xlabel('std mag',fontsize=12)
    ax.set_ylabel('residuals [mag]',fontsize=12)
    ax.plot(magstd, res, 'b.')
    ax.grid(True)
    plt.savefig(outputFile, format='png')
    
    outputFile = '%s.res_vs_color.png' % (outBaseName)
    if verbose > 1:  print 'Outputting '+outputFile
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('residuals vs std color')
    ax.set_xlabel('std color',fontsize=12)
    ax.set_ylabel('residuals [mag]',fontsize=12)
    ax.plot(colorstd, res, 'b.')
    ax.grid(True)
    plt.savefig(outputFile, format='png')
    
    outputFile = '%s.res_vs_ccd.png' % (outBaseName)
    if verbose > 1:  print 'Outputting '+outputFile
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('residuals vs ccd')
    ax.set_xlabel('ccd number',fontsize=12)
    ax.set_ylabel('residuals [mag]',fontsize=12)
    ax.plot(ccd, res, 'b.')
    ax.grid(True)
    plt.savefig(outputFile, format='png')
    
    outputFile = '%s.res_vs_expnum.png' % (outBaseName)
    if verbose > 1:  print 'Outputting '+outputFile
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('residuals vs expnum')
    ax.set_xlabel('exposure number',fontsize=12)
    ax.set_ylabel('residuals [mag]',fontsize=12)
    ax.plot(expnum, res, 'b.')
    ax.grid(True)
    plt.savefig(outputFile, format='png')
    
    outputFile = '%s.res_vs_mjd.png' % (outBaseName)
    if verbose > 1:  print 'Outputting '+outputFile
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('residuals vs mjd')
    ax.set_xlabel('MJD',fontsize=12)
    ax.set_ylabel('residuals [mag]',fontsize=12)
    ax.plot(mjdobs, res, 'b.')
    ax.grid(True)
    plt.savefig(outputFile, format='png')

    # Future:  add a section here to output a "residuals vs. X,Y position on the focal plane" plot...
    #outputFile = '%s.res_vs_fpXY.png' % (outBaseName)


    # Return with exit code 0...
    sys.exit(0)


#---------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#---------------------------------------------------------------------------
