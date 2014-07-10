#!/usr/bin/env python

# Authors:   Brian Yanny and Douglas Tucker
# Date:      17 May 2013
# Updated:   23 May 2013
# Updated:    9 Jul 2014

# Description:
#
# This file defines methods for solving the photometric zeropoints
# (a_1, a_2, ..., a_N), the instrumental color term coefficients
# (b_1, b_2, ..., b_N), and the first-order extinction (k) for a
# given band for a given night by fitting the following equation
# for standard star observations:
#    m_inst-m_std = a_1 + ... a_N +
#                   b_1*(stdColor-stdColor0) + ... +
#                   b_N*(stdColor-stdColor0) + kX,
# where m_inst is the instrumental (observed) mag of a standard
# star, m_std is the known calibrated mag of the standard star,
# stdColor is the known calibrated color of the standard star
# (e.g., its g-r color), stdColor0 is a zeropoint constant for
# the standard color and X is the airmass of the observation.
#
# For a camera with a single CCD, the above equation reduces to
# the following, simpler form:
#    m_inst-m_std = a + b*(stdColor-stdColor0) + kX
#
# For the explicit case of the g band and a g-r color,
# this single-CCD example looks like this:
#   g_inst-g_std = a + b*( (g-r) - (g-r)_0 ) + kX

import numpy
import sys
import math
import os
import getopt
import pyfits
import datetime

#---------------------------------------------------------------------------
# Client usage.

def usage():
    clientName = os.path.basename(sys.argv[0])
    print
    print 'Usage:'
    print ' %s <inmatches> <outak> <bandid> <niter> <thresholdit> [--ksolve] [--bsolve] [--verbose=0 (default)] [-h,--help]' % clientName
    print 'where:'
    print '   inmatches                       is the input match file                                 (required)'
    print '   outak                           is the output file                                      (required)'
    print '   bandid                          is the band id number (u=0,g=1,r=2,i=3,z=4,Y=5)         (required)'
    print '   niter                           is the number of iterations for the outlier rejection   (required)'
    print '   thresholdit                     is the threshold (in mag) of for the outlier rejection  (required)'
    print '   --ksolve                        is a toggle to solve for the k term coefficient         (optional)'
    print '   --bsolve                        is a toggle to solve for the b term coefficients        (optional)'
    #print '   --nite=\'20130221\'               is the nite of observation                              (optional)'
    print '   --project=\'OPS\'                 is the project id                                       (optional)'
    print '   --psmversion=\'pyPSM_v0.1\'       is the version of pyPSM being used                      (optional)'
    print '   --verbose                       is the verbosity level (default=0)                      (optional)'
    print '   -h,--help                       is a toggle to print out this usage guide               (optional)'
    print
    print 'Examples:'
    print ' %s matchemup.csv psm_solve_Y.txt 5 4 0.1 --bsolve --ksolve --verbose=3' % clientName
    print ' %s --help' % clientName
    print
    
    # value of nite is now read from the input match file
    #nite = '20130221'               # nite of observation
    project = 'OPS'                 # project name
    psmversion = 'pyPSM_v0.1'       # psm version id
    mag_type = 'mag_psf'            # type of magnitude used in fit (e.g., mag_psf, mag_aper_8, ...)

#---------------------------------------------------------------------------

def psm(inmatches,outak,bandid,niter,thresholdit,ksolve,bsolve,nite,project,psmversion,mag_type,verbose):
    
    if verbose>1:
        print 'psm arguments:'
        print '----------------------------------'
        print    'inmatches:      '+str(inmatches)
        print    'outak:          '+str(outak)
        print    'bandid:         '+str(bandid)
        print    'niter:          '+str(niter)
        print    'thresholdit:    '+str(thresholdit)
        print    'ksolve:         '+str(ksolve)
        print    'bsolve:         '+str(bsolve)
        #print    'nite:           '+str(nite)
        print    'project:        '+str(project)
        print    'psmversion:     '+str(psmversion)
        print    'mag_type:       '+str(mag_type)
        print    'verbose:        '+str(verbose)
        print
    #endif
    
    # Currently, pyPSM always solves for the a coefficients.
    # If, in the future, the option is added to choose whether
    # or not the a coefficients will be solved for, an "--asolve"
    # toggle should be added to Main method...
    asolve = 1
    
    # Use the signal value to denote values that are bad or meaningless...
    signalValue = -9999.00
    
    # Extract info from the list of arguments passed to psm.py...
    # First, convert certain arguments from strings to numbers...
    bandid = int(bandid)
    niter = int(niter)
    thresholdit = float(thresholdit)
    
    # Determine which band is to be solved for
    bandlist = ['u','g','r','i','z','Y']
    if bandid < len(bandlist):
        fitband = bandlist[bandid]
    else:
        print 'Bandid %d has no corresponding band...' % bandid
        sys.exit(1)
    #endif
    
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
    
    
    # The input matched file is inmatches,
    # the output matched file (used for iterative clipping) is inmatches.<fitband>.tmp,
    # the residuals file (used for QA plots) is inmatches.<fitband>.res, and
    # the FITS table (for upload to the NCSA DES db) is inmatches.<fitband>.fit.
    infile = inmatches
    outfile = inmatches+'.'+fitband+'.tmp'
    resfile = inmatches+'.'+fitband+'.res.csv'
    outfitsfile = inmatches+'.'+fitband+'.fit'
    usedcatfilesfile = inmatches+'.'+fitband+'.cats.list'
    
    ## The default zeropoint...
    #defaultzp = 0.0
    
    # The base magnitude error (for adding in quadrature with the Poisson error)
    baseMagErr = 0.02
    
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
        nite = -1
        
        # initialize list of image catalog files used in the PSM fit...
        usedcatfilesList = []
        
        # Open input file, read header, and identify necessary columns...
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
            if hn[i].upper() == 'MAGSTD':
                magstdcol=i
            if hn[i].upper() == 'COLORSTD':
                colorstdcol=i
        #endfor
        
        
        # Read input file line-by-line and populate AA and BB matrices
        # (to solve matrix equation AA.XX=BB)...
        for l in fd:
            lsp=l.strip().split(',')
            band=lsp[bandcol]
            if band != fitband:
                continue
            #endif
            
            colorstd=float(lsp[colorstdcol])
            if ((band == 'i' or band == 'z') and (colorstd < 0.0 or colorstd > 0.7)) or ((band == 'g' or band == 'r') and (colorstd < 0.2 or colorstd > 1.2)):
                continue
            #endif
            
            ccd = int(lsp[ccdcol])
            magstd = float(lsp[magstdcol])
            mag = float(lsp[magcol])
            magerr = float(lsp[magerrcol])
            #totMagErr = math.sqrt(baseMagErr*baseMagErr + magerr*mag)
            totMagErr = baseMagErr
            weight = 1.0/(baseMagErr*baseMagErr)
            exptime = float(lsp[exptimecol])
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
        ofd3=open(resfile,'w')
        ofd3.write('res,airmass,magstd,colorstd,ccd,expnum,mjdobs,crpix1,crpix2,naxis1,naxis2,ximage,yimage,catfilename\n')
        ofd2=open(outfile,'w')
        fd=open(infile)
        for l in fd:
            lsp=l.strip().split(',')
            band=lsp[bandcol]
            if band != fitband:
                continue
            
            colorstd=float(lsp[colorstdcol])
            if ((band == 'i' or band == 'z') and (colorstd < 0.0 or colorstd > 0.7)) or ((band == 'g' or band == 'r') and (colorstd < 0.2 or colorstd > 1.2)):
                continue
            
            ccd=int(lsp[ccdcol])
            magstd=float(lsp[magstdcol])
            mag=float(lsp[magcol])
            exptime=float(lsp[exptimecol])
            airmass=float(lsp[airmasscol])
            expnum=int(lsp[expnumcol])
            nite=lsp[nitecol]
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
        infile = outfile
        outfile = infile+'.tmp'
        
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
    
    # Find unique filenames in usedcatfilesList (using python "set" command),
    #  re-convert the set back to a list, sort the final and unique list, and
    #  output to (usedcatfilesfile...
    usedcatfilesSet = set(usedcatfilesList)
    usedcatfilesList = list (usedcatfilesSet)
    usedcatfilesList.sort()
    ofd4=open(usedcatfilesfile,'w')
    ofd4.write('FILENAME\n')
    for i in range(0,len(usedcatfilesList)):
        ofd4.write(usedcatfilesList[i]+'\n')
    #endfor
    ofd4.close()
    
    
    # Output the results of fit...
    os.system("/bin/rm -f "+inmatches+'.'+fitband+'.tmp*')
    ofd=open(outak,'w')
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
    os.system("/bin/rm -f "+outfitsfile)
    ccdidArray = numpy.arange(1,63)
    niteFormat = 'A%d' % len(nite)
    niteArray = numpy.array([nite]*62)
    mjdloArray = mjdobsLo*numpy.ones(nccd,dtype=numpy.int)
    mjdhiArray = mjdobsHi*numpy.ones(nccd,dtype=numpy.int)
    ccdidArray = numpy.arange(1,63)
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
    col4  = pyfits.Column(name='ccdid',          format='I', array=ccdidArray)
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
    
    
    print "That's all, folks!"
    
    sys.exit(0)

#---------------------------------------------------------------------------
# Main method:

if __name__ == "__main__":
    
    # If insufficient number of required arguments, print out usage...
    if len(sys.argv[1:]) < 5:
        usage()
        sys.exit(1)
    #endif
    
    # Parse required arguments (arguments 1-5)...
    inmatches   = sys.argv[1]
    outak       = sys.argv[2]
    bandid      = sys.argv[3]
    niter       = sys.argv[4]
    thresholdit = sys.argv[5]
    
    # Parse any optional arguments (any beyond argument 5)...
    try:
        opts,args = getopt.getopt(sys.argv[6:],'hv',['bsolve', 'ksolve', 'nite=', 'project=', 'psmversion=', 'mag_type=', 'help', 'verbose='])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    #end try
    
    # Default values for the optional arguments...
    verbose = 0                     # Verbosity level (amount of output to screen): 0=little; >3=lots
    bsolve = 0                      # Solve for b coefficients? (0=no, 1=yes)
    ksolve = 0                      # Solve for k coefficient?  (0=no, 1=yes)
    nite = '20130221'               # nite of observation
    project = 'OPS'                 # project name
    psmversion = 'pyPSM_v0.1'       # psm version id
    mag_type = 'mag_psf'            # type of magnitude used in fit (e.g., mag_psf, mag_aper_8, ...)
    
    for o, a in opts:
        if o == '--bsolve':
            bsolve = 1
        elif o == '--ksolve':
            ksolve = 1
        elif o == '--nite':
            nite = a
        elif o == '--project':
            project = a
        elif o == '--psmversion':
            psmversion = a
        elif o == '--mag_type':
            mag_type = a
        elif o in ('-h', '--help'):
            usage()
            sys.exit(0)
        elif o in ('-v', '--verbose'):
            verbose = int(a)
    #endif
    #endfor
    
    # Call psm method
    psm(inmatches,outak,bandid,niter,thresholdit,ksolve,bsolve,nite,project,psmversion,mag_type,verbose)


#---------------------------------------------------------------------------


