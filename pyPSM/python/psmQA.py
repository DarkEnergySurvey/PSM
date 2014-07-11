#!/usr/bin/env python

# Author:    Douglas Tucker
# Date:      21 May 2013
# Updated:   21 May 2013

# Description:
#
# This file defines methods for creating QA plots from the residuals
# residuals CSV file output by psm.

import numpy
import sys
import math
import os
import getopt
import csv
import matplotlib.pyplot as plt

#---------------------------------------------------------------------------
# Client usage.

def usage():
    clientName = os.path.basename(sys.argv[0])
    print
    print 'Usage:'
    print ' %s <inmatches> <outak> <bandid> <niter> <thresholdit> [--ksolve] [--bsolve] [--verbose=0 (default)] [-h,--help]' % clientName
    print 'where:'
    print '   psmQATableFile is the input residuals file (in CSV format)          (required)'
    print '   --verbose    is the verbosity level (default=0)                     (optional)'
    print '   -h,--help    is a toggle to print out this usage guide              (optional)'
    print
    print 'Examples:'
    print ' %s matchemup.g.res.csv --verbose=3' % clientName
    print ' %s --help' % clientName
    print


#---------------------------------------------------------------------------

def psmQA(psmQATableFile):
    
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
    
    # Output QA plots...
    
    outputFile = 'PSM_QA_res_vs_airmass.png'
    print 'Outputting '+outputFile
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('residuals vs airmass')
    ax.set_xlabel('airmass',fontsize=12)
    ax.set_ylabel('residuals [mag]',fontsize=12)
    ax.plot(airmass, res, 'b.')
    ax.grid(True)
    plt.savefig(outputFile, format='png')
    
    outputFile = 'PSM_QA_res_vs_mag.png'
    print 'Outputting '+outputFile
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('residuals vs std mag')
    ax.set_xlabel('std mag',fontsize=12)
    ax.set_ylabel('residuals [mag]',fontsize=12)
    ax.plot(magstd, res, 'b.')
    ax.grid(True)
    plt.savefig(outputFile, format='png')
    
    outputFile = 'PSM_QA_res_vs_color.png'
    print 'Outputting '+outputFile
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('residuals vs std color')
    ax.set_xlabel('std color',fontsize=12)
    ax.set_ylabel('residuals [mag]',fontsize=12)
    ax.plot(colorstd, res, 'b.')
    ax.grid(True)
    plt.savefig(outputFile, format='png')
    
    outputFile = 'PSM_QA_res_vs_ccd.png'
    print 'Outputting '+outputFile
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('residuals vs ccd')
    ax.set_xlabel('ccd number',fontsize=12)
    ax.set_ylabel('residuals [mag]',fontsize=12)
    ax.plot(ccd, res, 'b.')
    ax.grid(True)
    plt.savefig(outputFile, format='png')
    
    outputFile = 'PSM_QA_res_vs_expnum.png'
    print 'Outputting '+outputFile
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('residuals vs expnum')
    ax.set_xlabel('exposure number',fontsize=12)
    ax.set_ylabel('residuals [mag]',fontsize=12)
    ax.plot(expnum, res, 'b.')
    ax.grid(True)
    plt.savefig(outputFile, format='png')
    
    
    
    print "That's all, folks!"
    
    sys.exit(0)

#---------------------------------------------------------------------------
# Main method:

if __name__ == "__main__":
    
    # If help requested, or if insufficient number of arguments, print out usage.
    if (sys.argv[1] == '-h') or (sys.argv[1] == '--help'):
        usage()
        sys.exit(0)
    elif len(sys.argv[1:]) < 1:
        usage()
        sys.exit(1)
    #endif
    
    # Parse required argument (arguments 1)...
    psmQATableFile = sys.argv[1]
    
    # Parse any optional arguments (any beyond argument 1)...
    try:
        opts,args = getopt.getopt(sys.argv[2:],'h',['help', 'verbose='])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    #end try
    
    verbose = 0        # Default value for verbosity
    
    for o, a in opts:
        if o in ('-h', '--help'):
            usage()
            sys.exit(0)
        elif o in ('--verbose'):
            verbose = int(a)
    #endif
    #endfor
    
    # Call psmQA method
    psmQA(psmQATableFile)


#---------------------------------------------------------------------------


