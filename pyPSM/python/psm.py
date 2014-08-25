#!/usr/bin/env python

# Authors:   Brian Yanny and Douglas Tucker
# Date:      17 May 2013
# Updated:   23 May 2013
# Updated:   10 Jul 2014
# Updated:   18 Aug 2014
# Updated:   25 Aug 2014

"""
psm.py

Description:

This file defines methods for solving the photometric zeropoints
(a_1, a_2, ..., a_N), the instrumental color term coefficients
(b_1, b_2, ..., b_N), and the first-order extinction (k) for a
given band for a given night by fitting the following equation
for standard star observations:
   m_inst-m_std = a_1 + ... a_N +
                  b_1*(stdColor-stdColor0) + ... +
                  b_N*(stdColor-stdColor0) + kX,
where m_inst is the instrumental (observed) mag of a standard
star, m_std is the known calibrated mag of the standard star,
stdColor is the known calibrated color of the standard star
(e.g., its g-r color), stdColor0 is a zeropoint constant for
the standard color and X is the airmass of the observation.

For a camera with a single CCD, the above equation reduces to
the following, simpler form:
   m_inst-m_std = a + b*(stdColor-stdColor0) + kX

For the explicit case of the g band and a g-r color,
this single-CCD example looks like this:
   g_inst-g_std = a + b*( (g-r) - (g-r)_0 ) + kX


Examples:

psm.py --help

psm.py --inputMatchFile matched-20131002-g-r03p01.csv --outputResultsFITSFile psmResults-20131002-g-r03p01.fits --outputResultsLogFile psmResults-20131002-g-r03p01.log --outputResidualsFile psmResiduals-20131002-g-r03p01.csv --outputCatsUsedFile psmCats-20131002-g-r03p01.list --band g --niter 3 --thresholdit 0.1 --ksolve --bsolve --ccdExcludeList 61,2,31 --expnumExcludeList 240688,240724,240774 --verbose 3

"""

import numpy
import sys
import math
import os
import getopt
import pyfits
import datetime

#---------------------------------------------------------------------------

def main():
    import argparse
    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputMatchFile', help='CSV file containing the matches from matchsort.py', default='matched.csv')
    parser.add_argument('--outputResultsFITSFile', help='FITS table containing the results from the fit (for upload to the DESDM database)', default='psmResults.fits')
    parser.add_argument('--outputResultsLogFile', help='ASCII text file containing a human-readable log of the results from the fit', default='psmResults.log')
    parser.add_argument('--outputResidualsFile', help='CSV table containing the final fit residuals vs. a variety of parameters (used downstreadm by psmQA.py)', default='psmResiduals.csv')
    parser.add_argument('--outputCatsUsedFile', help='ASCII list of the image catalogs used in the final fit)', default='psmCats.list')
    parser.add_argument('--niter', help='the number of iterations for the outlier rejection', type=int, default=0)
    parser.add_argument('--thresholdit', help='the threshold (in mag) for the outlier rejection', type=float, default=0)
    parser.add_argument('--band', help='band to be fit', default='g', choices=['u','g','r','i','z','Y'])
    parser.add_argument('--bsolve',help='a toggle to solve for the b term coefficients', default=False, action='store_true')
    parser.add_argument('--ksolve',help='a toggle to solve for the k term coefficient', default=False, action='store_true')
    parser.add_argument('--expnumExcludeList', help='a comma-separated list of EXPNUMs to exclude from the fit', default='')
    parser.add_argument('--ccdExcludeList', help='a comma-separated list of CCDNUMs to exclude from the fit', default='')
    parser.add_argument('--exptimeLo', help='exposure time lower limit (in seconds) to use in fit', type=float, default=2.0)
    parser.add_argument('--exptimeHi', help='exposure time upper limit (in seconds) to use in fit', type=float, default=100.0)
    parser.add_argument('--magLo', help='std star mag lower limit  to use in fit', type=float, default=15.0)
    parser.add_argument('--magHi', help='std star mag upper limit to use in fit', type=float, default=18.0)
    parser.add_argument('--project', help='the project id (deprecated?)', default='OPS')
    parser.add_argument('--psmversion', help='the version of pyPSM being used (deprecated?)', default='pyPSM_v0.1')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0, 1, 2, ...)', type=int, default=0)

    args = parser.parse_args()

    if args.verbose > 0: print args

    psm(args)

#---------------------------------------------------------------------------

def psm(args):

    infile = args.inputMatchFile
    outfitsfile = args.outputResultsFITSFile
    outlogfile = args.outputResultsLogFile
    outresfile = args.outputResidualsFile
    outcatsfile = args.outputCatsUsedFile
    fitband = args.band
    niter = args.niter
    thresholdit = args.thresholdit
    exptimelo = args.exptimeLo
    exptimehi = args.exptimeHi
    maglo = args.magLo
    maghi = args.magHi
    project = args.project
    psmversion = args.psmversion
    verbose = args.verbose

    ccdExcludeArray = numpy.sort(numpy.fromstring(args.ccdExcludeList,dtype=int,sep=','))
    expnumExcludeArray = numpy.sort(numpy.fromstring(args.expnumExcludeList,dtype=int,sep=','))

    if args.ksolve:
        ksolve=1
    else:
        ksolve=0

    if args.bsolve:
        bsolve=1
    else:
        bsolve=0    
    
    # Currently, pyPSM always solves for the a coefficients.
    # If, in the future, the option is added to choose whether
    # or not the a coefficients will be solved for, an "--asolve"
    # toggle should be added to Main method...
    asolve = 1

    if verbose>1:
        print 'psm.py parameters:'
        print '----------------------------------'
        print    'infile:           '+str(infile)
        print    'outfitsfile:      '+str(outfitsfile)
        print    'outlogfile:       '+str(outlogfile)
        print    'outresfile:       '+str(outresfile)
        print    'outcatsfile:      '+str(outcatsfile)
        print    'fitband:          '+str(fitband)
        print    'niter:            '+str(niter)
        print    'thresholdit:      '+str(thresholdit)
        print    'asolve:           '+str(asolve)
        print    'bsolve:           '+str(bsolve)
        print    'ksolve:           '+str(ksolve)
        print    'ccdExcludeList    '+str(ccdExcludeArray)
        print    'expnumExcludeList '+str(expnumExcludeArray)
        print    'exptimelo         '+str(exptimelo)
        print    'exptimehi         '+str(exptimehi)
        print    'maglo             '+str(maglo)
        print    'maghi             '+str(maghi)
        print    'project:          '+str(project)
        print    'psmversion:       '+str(psmversion)
        print    'verbose:          '+str(verbose)
        print
    #endif
    
    # Use the signal value to denote values that are bad or meaningless...
    signalValue = -9999.00
    
    # Set certain defaults for each band...
    if fitband == 'u':
        cfitband = 'g'
        stdColorName = 'u-g'
        bdefault = 0.06
        color0 = 1.39
        kdefault = 0.489
        # The u-band kgoodLo and kgoodHi need verification...
        kgoodLo = 0.2
        kgoodHi = 0.8
    elif fitband == 'g':
        cfitband = 'r'
        stdColorName = 'g-r'
        bdefault = -0.1
        color0 = 0.53
        kdefault = 0.181
        kgoodLo = 0.15
        kgoodHi = 0.30
    elif fitband == 'r':
        cfitband = 'g'
        stdColorName = 'g-r'
        bdefault = -0.08
        color0 = 0.53
        kdefault = 0.095
        kgoodLo = 0.06
        kgoodHi = 0.20
    elif fitband == 'i':
        cfitband = 'z'
        stdColorName = 'i-z'
        bdefault = -0.3
        color0 = 0.09
        kdefault = 0.089
        kgoodLo = 0.03
        kgoodHi = 0.16
    elif fitband == 'z':
        cfitband = 'i'
        stdColorName = 'i-z'
        bdefault = -0.09
        color0 = 0.053
        kdefault = 0.089
        kgoodLo = 0.03
        kgoodHi = 0.15
    elif fitband == 'Y':
        cfitband = 'z'
        stdColorName = 'z-Y'
        bdefault = 0.26
        color0 = 0.05
        kdefault = 0.050
        kgoodLo = 0.01
        kgoodHi = 0.17
    #endif
    
    ## The default zeropoint...
    #defaultzp = 0.0
    
    # The base magnitude error (for adding in quadrature with the Poisson error)
    baseMagErr = 0.02

    # Temporary file, to contain culled matched data after each iteration...
    tmpfile = infile+'.tmp'
    
    # The number of parameters as a whole (nparam), and
    # the number of free parameters (nFreeParam), which
    # depends on whether or not we are solving for the
    # b and k coefficients...
    nccd = 62
    nparam = 2*nccd + 1
    nFreeParam = nccd
    if ksolve==1:
        nFreeParam = nFreeParam + 1
    if bsolve==1:
        nFreeParam = nFreeParam + nccd
    if verbose > 1:
        print ksolve, bsolve, nparam, nFreeParam
    
    # Iteration loop...
    for iiter in range(0,niter):
        # Initialize various variables used in the fit,
        # including the matrices AA and BB (for the
        # matrix equation AA.XX=BB)...
        sum = 0.0
        sumsq = 0.0
        sumchi = 0.0
        sumchisq = 0.0
        ninfit = 0.0
        #zv = numpy.zeros(nparam)
        II = numpy.array(numpy.identity(nparam))
        for i in range(0,nparam):
            for j in range(0,nparam):
                II[i][j] = 0
        BB = numpy.zeros(nparam)
        #AA = numpy.vstack([zv,II]).T
        AA = II
        
        # counter of stars on each CCD...
        nstar = numpy.zeros(nccd,dtype=numpy.int)
        
        # initialize range of mjd's covered by the standard star observations included in the fit
        mjdobsLo =  1.00e6
        mjdobsHi = -1.00e6
        
        # initialize value of nite; we currently assume that there is only one nite represented
        #   in the input match file...
        # (should add sanity check at some point to verify there is only only nite represented)
        nite = ''

        # initialize value of mag_type; we currently assume that there is only one mag_type represented
        #   in the input match file...
        # (should add sanity check at some point to verify there is only only mag_type represented)
        mag_type = ''

        # initialize list of image catalog files used in the PSM fit...
        usedcatfilesList = []
        
        # Open input file, read header, and identify necessary columns...
        # (should be able to move this subsection to outside the iteration loop)
        fd=open(infile)
        h=fd.readline()
        hn=h.strip().split(',')
        for i in range(0,len(hn)):
            if hn[i].upper() == 'NITE':
                nitecol=i
            if hn[i].upper() == 'MJD_OBS':
                mjdobscol=i
            if hn[i].upper() == 'EXPNUM':
                expnumcol=i
            if hn[i].upper() == 'BAND':
                bandcol=i
            if hn[i].upper() == 'EXPTIME':
                exptimecol=i
            if hn[i].upper() == 'AIRMASS':
                airmasscol=i
            if hn[i].upper() == 'CCDNUM':
                ccdcol=i
            if hn[i].upper() == 'FILENAME':
                catfilenamecol=i
            if hn[i].upper() == 'CRPIX1':
                crpix1col=i
            if hn[i].upper() == 'CRPIX2':
                crpix2col=i
            if hn[i].upper() == 'NAXIS1':
                naxis1col=i
            if hn[i].upper() == 'NAXIS2':
                naxis2col=i
            if hn[i].upper() == 'X_IMAGE':
                ximagecol=i
            if hn[i].upper() == 'Y_IMAGE':
                yimagecol=i
            if hn[i].upper() == 'MAG':
                magcol=i
            if hn[i].upper() == 'MAGERR':
                magerrcol=i
            if hn[i].upper() == 'ZEROPOINT':
                zeropointcol=i
            if ( (hn[i].upper() == 'MAGTYPE') or (hn[i].upper() == 'MAG_TYPE') ):
                magtypecol=i
            if hn[i].upper() == 'MAGSTD':
                magstdcol=i
            if hn[i].upper() == 'COLORSTD':
                colorstdcol=i
        #endfor
        
        
        # Read input file line-by-line and populate AA and BB matrices
        # (to solve matrix equation AA.XX=BB)...
        for l in fd:

            lsp=l.strip().split(',')

            # skip if the wrong filter band...
            band=lsp[bandcol]
            if band != fitband:
                continue
            #endif
            
            exptime = float(lsp[exptimecol])
            if ( (exptime < exptimelo) or (exptime > exptimehi) ):
                continue
            #endif

            # skip if on the expnum exclude list...
            expnum=int(lsp[expnumcol])
            if expnum in expnumExcludeArray:
                continue
            #endif

            # skip if on the ccd exclude list...
            ccd = int(lsp[ccdcol])
            if ccd in ccdExcludeArray:
                continue
            #endif

            # skip if outside the preferred std star mag range...
            magstd = float(lsp[magstdcol])
            if ( (magstd < maglo) or (magstd > maghi) ):
                continue
            #endif
            
            # skip if outside the preferred std star color range...
            colorstd=float(lsp[colorstdcol])
            if ((band == 'i' or band == 'z') and (colorstd < 0.0 or colorstd > 0.7)) or ((band == 'g' or band == 'r') and (colorstd < 0.2 or colorstd > 1.2)):
                continue
            #endif
            
            mag = float(lsp[magcol])
            magerr = float(lsp[magerrcol])
            #totMagErr = math.sqrt(baseMagErr*baseMagErr + magerr*mag)
            totMagErr = baseMagErr
            weight = 1.0/(baseMagErr*baseMagErr)
            airmass = float(lsp[airmasscol])
            zeropoint = float(lsp[zeropointcol])
            dm  = (mag - (zeropoint-25.) + 2.5*math.log(exptime,10)) - magstd
            #dm  = (mag - defaultzp + 2.5*math.log(exptime,10)) - magstd
            
            if abs(dm) < 4.5:
                
                nstar[ccd-1] += 1
                
                #Code to indices:
                # k:       index 0
                # a[ccd]:  index iparam_a
                # b[ccd]:  index iparam_b
                
                iparam_a = ccd
                iparam_b = nccd + ccd
                
                #Slight difference in weighting between original pyPSM version and javaPSM version.
                #Need to check which is correct.
                #Original pyPSM version:
                #BB[0] += weight*weight*airmass*dm
                #BB[iparam_a] += weight*dm
                #AA[0][0] += weight*weight*airmass*airmass
                #AA[iparam_a][0] += weight*airmass
                #AA[0][iparam_a] += weight*airmass
                #AA[iparam_a][iparam_a] += 1.0*weight
                
                #Version from javaPSM (note BB[0] and A[0][0] have each lost one factor of weight):
                BB[0] += weight*airmass*dm
                BB[iparam_a] += weight*dm
                AA[0][0] += weight*airmass*airmass
                AA[iparam_a][0] += weight*airmass
                AA[0][iparam_a] += weight*airmass
                AA[iparam_a][iparam_a] += 1.0*weight
                
                dc = colorstd-color0
                
                if bsolve==1:
                    AA[iparam_b][iparam_b] += dc * dc * weight
                    AA[0][iparam_b] += dc * airmass * weight
                    AA[iparam_b][0] += dc * airmass * weight
                    AA[iparam_a][iparam_b] += dc * weight
                    AA[iparam_b][iparam_a] += dc * weight
                    BB[iparam_b] += dc * dm * weight
                else:
                    AA[iparam_b][iparam_b] = 1.0
                    AA[0][iparam_b] += dc * airmass * weight
                    AA[iparam_b][0] = 0.0
                    AA[iparam_a][iparam_b] += dc * weight
                    AA[iparam_b][iparam_a] = 0.0
                    #BB[iparam_b] = bdefaultValues[iccd]
                    BB[iparam_b] = bdefault
                #endif
        
            #endif
        
        #endfor
        
        # Find any CCDs for which nstar==0
        badinds = numpy.hstack(numpy.nonzero(nstar==0))
        
        # Set the a's, b's of any CCDs for which nstar=0 to signal values,
        # reducing the number of free parameters as appropriate
        for iccd in range(0,len(badinds)):
            
            ccd = badinds[iccd] + 1
            iparam_a = ccd
            iparam_b = nccd + ccd
            
            # Fix a's...
            AA[iparam_a,:] = 0.0
            AA[iparam_a][iparam_a] = 1.0
            BB[iparam_a] = signalValue+25.00
            # Since we always solve for a (asolve always equals 1),
            # any CCD with nstar=0 reduces the number of free
            # parameters by 1...
            if asolve == 1:
                nFreeParam += -1
            #endif
            
            # Fix b's...
            AA[iparam_b,:] = 0.0
            AA[iparam_b][iparam_b] = 1.0
            BB[iparam_b] = signalValue
            # Only reduce the number of free parameters if we are
            # solving for the b term coefficients...
            if bsolve == 1:
                nFreeParam += -1
            #endif
        
        #endfor
        
        if (ksolve==0):
            AA[0][0] = 1.0
            BB[0] = kdefault
            for iccd in range(0,int(nccd)):
                ccd = iccd+1
                iparam_a = ccd
                iparam_b = nccd + ccd
                AA[0][iparam_a] = 0.0
                AA[0][iparam_b] = 0.0
            #endfor
        #endif
        
        # Solve for XX in AA.XX=BB matrix equation...
        #(XX,residue,rank,s) = numpy.linalg.lstsq(AA,BB)
        #print 'rank', rank
        AAinv = numpy.linalg.inv(AA)
        XX=AAinv.dot(BB)
        # Note:  verify that these errors are what we expect:
        errors = numpy.sqrt(numpy.diagonal(AAinv))
        # Reset aerr and berr to the signal value for any CCDs for which nstar=0...
        for iccd in range(0,len(badinds)):
            ccd = badinds[iccd] + 1
            iparam_a = ccd
            iparam_b = nccd + ccd
            errors[iparam_a] = signalValue
            errors[iparam_b] = signalValue
        #endfor
        
        # Output to matched catalog output file those stars that
        # survive this iteration's clipping...
        ofd3=open(outresfile,'w')
        ofd3.write('res,airmass,magstd,colorstd,ccd,expnum,mjdobs,crpix1,crpix2,naxis1,naxis2,ximage,yimage,catfilename\n')
        ofd2=open(tmpfile,'w')
        fd=open(infile)
        for l in fd:

            lsp=l.strip().split(',')

            # skip if the wrong filter band...
            band=lsp[bandcol]
            if band != fitband:
                continue
            #endif
            
            exptime = float(lsp[exptimecol])
            if ( (exptime < exptimelo) or (exptime > exptimehi) ):
                continue
            #endif

            # skip if on the expnum exclude list...
            expnum=int(lsp[expnumcol])
            if expnum in expnumExcludeArray:
                continue
            #endif

            # skip if on the ccd exclude list...
            ccd = int(lsp[ccdcol])
            if ccd in ccdExcludeArray:
                continue
            #endif

            # skip if outside the preferred std star mag range...
            magstd = float(lsp[magstdcol])
            if ( (magstd < maglo) or (magstd > maghi) ):
                continue
            #endif
            
            # skip if outside the preferred std star color range...
            colorstd=float(lsp[colorstdcol])
            if ((band == 'i' or band == 'z') and (colorstd < 0.0 or colorstd > 0.7)) or ((band == 'g' or band == 'r') and (colorstd < 0.2 or colorstd > 1.2)):
                continue
            #endif
            
            mag=float(lsp[magcol])
            airmass=float(lsp[airmasscol])
            nite=lsp[nitecol]
            mag_type=lsp[magtypecol]
            mjdobs=float(lsp[mjdobscol])
            crpix1=float(lsp[crpix1col])
            crpix2=float(lsp[crpix2col])
            naxis1=int(lsp[naxis1col])
            naxis2=int(lsp[naxis2col])
            ximage=float(lsp[ximagecol])
            yimage=float(lsp[yimagecol])
            catfilename=lsp[catfilenamecol].strip()
            # Note:  should change befault to XX[nccd+ccd]
            #dm = (mag + 2.5*math.log(exptime,10) - XX[ccd] - XX[0]*airmass - bdefault*(colorstd-color0)) - magstd
            dm = (mag + 2.5*math.log(exptime,10) - XX[ccd] - XX[0]*airmass - XX[nccd+ccd]*(colorstd-color0)) - magstd
            sum += dm
            sumsq += dm*dm
            sumchi += dm/totMagErr
            sumchisq += (dm/totMagErr)*(dm/totMagErr)
            ninfit += 1
            # Should change to a sigma-clip...
            if abs(dm) < thresholdit:
                if (mjdobs < mjdobsLo):  mjdobsLo = mjdobs
                if (mjdobs > mjdobsHi):  mjdobsHi = mjdobs
                usedcatfilesList.append(catfilename)
                ofd2.write(l)
                resOutputLine = '%.4f,%.3f,%.4f,%.4f,%d,%d,%.5f,%.3f,%.3f,%d,%d,%.5f,%.5f,%s\n' % (dm,airmass,magstd,colorstd,ccd,expnum,mjdobs,crpix1,crpix2,naxis1,naxis2,ximage,yimage,catfilename)
                ofd3.write(resOutputLine)
            #endif
        
        #endfor (fd)
        
        ofd3.close()
        ofd2.close()
        infile = tmpfile
        tmpfile = infile+'.tmp'
        
        # Calculate characteristics of the fit...
        photometricFlag = -1
        dof = ninfit - nFreeParam
        if dof > 0:
            avg = float(sum/ninfit)
            sigma = math.sqrt(float(ninfit)/float(dof))*math.sqrt(sumsq/ninfit-avg*avg)
            avgchi = float(sumchi/ninfit)
            chisq = math.sqrt(float(ninfit)/float(dof))*math.sqrt(sumchisq/ninfit-avgchi*avgchi)
            if (sigma < 0.025) and (XX[0] >= kgoodLo and XX[0] <= kgoodHi):
                photometricFlag = 1
            else:
                photometricFlag = 0
            #endif
            if verbose > 3:
                print 'avg:'+str(avg)+'  avgchi:'+str(avgchi)+'  chisq:'+str(chisq)
                print 'k:'+str(XX[0])+' rms:'+str(sigma)+' n:'+str(ninfit)
                print 'kerr:'+str(errors[0])
            #endif
        else:
            print 'Number of free parameters ('+str(nFreeParam)+') >= number of data points ('+str(ninfit)+')'
            print 'Exiting now!'
            sys.exit(1)
        #endif
    
    #endfor (iiter)

    # Delete old tmp files...
    tmpfile = args.inputMatchFile
    for iiter in range(0,niter):
        tmpfile = tmpfile+'.tmp'
        os.system('/bin/rm -f '+tmpfile)
    #endfor


    # Find unique filenames in usedcatfilesList (using python "set" command),
    #  re-convert the set back to a list, sort the final and unique list, and
    #  output to (outcatsfile...
    usedcatfilesSet = set(usedcatfilesList)
    usedcatfilesList = list (usedcatfilesSet)
    usedcatfilesList.sort()
    ofd4=open(outcatsfile,'w')
    #removing to be compatible with workflow list no-header convention
    #ofd4.write('FILENAME\n')
    for i in range(0,len(usedcatfilesList)):
        ofd4.write(usedcatfilesList[i]+'\n')
    #endfor
    ofd4.close()
    
    
    # Output the results of fit...
    ofd=open(outlogfile,'w')
    ofd.write('niter_'+band+' '+str('%d'%int(niter))+'\n')
    ofd.write('thresholdit_'+band+' '+str('%.4f'%float(thresholdit))+'\n')
    ofd.write('n_'+band+' '+str('%d'%int(ninfit))+'\n')
    ofd.write('dof_'+band+' '+str('%d'%int(dof))+'\n')
    ofd.write('rms_'+band+' '+str('%.4f'%float(sigma))+'\n')
    ofd.write('chisq_'+band+' '+str('%.4f'%float(chisq))+'\n')
    ofd.write('k_'+band+' '+str('%.3f'%float(XX[0]))+' +/- '+str('%.3f'%float(errors[0]))+'\n')
    for i in range(1,63):
        outputLine = 'a_%s %2d %.3f +/- %.3f \n' % (band, i, XX[i]-25., errors[i])
        ofd.write(outputLine)
        if verbose > 1: print outputLine,
    #endfor
    for i in range(1,63):
        outputLine = 'b_%s %2d %.3f +/- %.3f \n' % (band, i, XX[nccd+i], errors[nccd+i])
        ofd.write(outputLine)
        if verbose > 1: print outputLine,
    #endfor
    ofd.close()
    
    # Output the results of fit into a FITS table file...
    os.system('/bin/rm -f '+outfitsfile)
    ccdnumArray = numpy.arange(1,63)
    niteFormat = 'A%d' % len(nite)
    niteArray = numpy.array([nite]*62)
    mjdloArray = mjdobsLo*numpy.ones(nccd,dtype=numpy.int)
    mjdhiArray = mjdobsHi*numpy.ones(nccd,dtype=numpy.int)
    ccdnumArray = numpy.arange(1,63)
    bandFormat = 'A%d' % len(band)
    bandArray = numpy.array([band]*62)
    kArray = XX[0]*numpy.ones(nccd,dtype=numpy.int)
    kerrArray = errors[0]*numpy.ones(nccd,dtype=numpy.int)
    aArray = XX[1:nccd+1]-25.
    aerrArray = errors[1:nccd+1]
    bArray = XX[nccd+1:]
    berrArray = errors[nccd+1:]
    rmsArray = sigma*numpy.ones(nccd,dtype=numpy.int)
    chisqArray = chisq*numpy.ones(nccd,dtype=numpy.int)
    dofArray = dof*numpy.ones(nccd,dtype=numpy.int)
    photometricFlagArray = photometricFlag*numpy.ones(nccd,dtype=numpy.int)
    psmversionFormat = 'A%d' % len(psmversion)
    psmversionArray = numpy.array([psmversion]*62)
    fit_timestamp = str(datetime.datetime.utcnow())
    fit_timestampFormat = 'A%d' % len(fit_timestamp)
    fit_timestampArray = numpy.array([fit_timestamp]*62)
    cbandFormat = 'A%d' % len(cfitband)
    cbandArray = numpy.array([cfitband]*62)
    stdcolor0Array = color0*numpy.ones(nccd,dtype=numpy.int)
    asolveArray = asolve*numpy.ones(nccd,dtype=numpy.int)
    bsolveArray = bsolve*numpy.ones(nccd,dtype=numpy.int)
    ksolveArray = ksolve*numpy.ones(nccd,dtype=numpy.int)
    projectFormat = 'A%d' % len(project)
    projectArray = numpy.array([project]*62)
    mag_typeFormat = 'A%d' % len(mag_type)
    mag_typeArray = numpy.array([mag_type]*62)
    
    col1  = pyfits.Column(name='nite',           format=niteFormat, array=niteArray)
    col2  = pyfits.Column(name='mjdlo',          format='D', array=mjdloArray)
    col3  = pyfits.Column(name='mjdhi',          format='D', array=mjdhiArray)
    col4  = pyfits.Column(name='ccdnum',         format='I', array=ccdnumArray)
    col5  = pyfits.Column(name='band',           format=bandFormat, array=bandArray)
    col6  = pyfits.Column(name='a',              format='E', array=aArray)
    col7  = pyfits.Column(name='aerr',           format='E', array=aerrArray)
    col8  = pyfits.Column(name='b',              format='E', array=bArray)
    col9 = pyfits.Column(name='berr',           format='E', array=berrArray)
    col10 = pyfits.Column(name='k',              format='E', array=kArray)
    col11 = pyfits.Column(name='kerr',           format='E', array=kerrArray)
    col12 = pyfits.Column(name='rms',            format='E', array=rmsArray)
    col13 = pyfits.Column(name='chi2',           format='E', array=chisqArray)
    col14 = pyfits.Column(name='dof',            format='I', array=dofArray)
    col15 = pyfits.Column(name='photomtricflag', format='I', array=photometricFlagArray)
    col16 = pyfits.Column(name='psmversion',     format=psmversionFormat, array=psmversionArray)
    col17 = pyfits.Column(name='fit_timestamp',  format=fit_timestampFormat, array=fit_timestampArray)
    col18 = pyfits.Column(name='cband',          format=cbandFormat, array=cbandArray)
    col19 = pyfits.Column(name='stdcolor0',      format='E', array=stdcolor0Array)
    col20 = pyfits.Column(name='asolve',         format='I', array=asolveArray)
    col21 = pyfits.Column(name='bsolve',         format='I', array=bsolveArray)
    col22 = pyfits.Column(name='ksolve',         format='I', array=ksolveArray)
    col23 = pyfits.Column(name='project',        format=projectFormat, array=projectArray)
    col24 = pyfits.Column(name='mag_type',       format=mag_typeFormat, array=mag_typeArray)
    
    cols=pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21, col22, col23, col24])
    tbhdu=pyfits.new_table(cols)
    tbhdu.writeto(outfitsfile)
    
    sys.exit(0)

#---------------------------------------------------------------------------

if __name__ == "__main__":
    main()

#---------------------------------------------------------------------------
