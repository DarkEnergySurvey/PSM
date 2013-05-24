#!/usr/bin/env python

# Authors:   Brian Yanny and Douglas Tucker
# Date:      17 May 2013
# Updated:   23 May 2013

# Description:
#
# This file defines methods for solving the photometric zeropoints 
# (a_1, a_2, ..., a_N), the instrumental color term coefficients
# (b_1, b_2, ..., b_N), and the first-order extinction (k) for a 
# given filter for a given night by fitting the following equation 
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
# For the explicit case of the g filter and a g-r color,
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
   print '   inmatches    is the input match file                                (required)'
   print '   outak        is the output file                                     (required)'
   print '   bandid       is the band id number (u=0,g=1,r=2,i=3,z=4,Y=5)        (required)'
   print '   niter        is the number of iterations for the outlier rejection  (required)'
   print '   thresholdit  is the threshold (in mag) of for the outlier rejection (required)'
   print '   --ksolve     is a toggle to solve for the k term coefficient        (optional)'
   print '   --ksolve     is a toggle to solve for the b term coefficients       (optional)'
   print '   --verbose    is the verbosity level (default=0)                     (optional)'
   print '   -h,--help    is a toggle to print out this usage guide              (optional)'
   print 
   print 'Examples:'
   print ' %s matchemup.csv psm_solve_Y.txt 5 4 0.1 --bsolve --ksolve --verbose=3' % clientName
   print ' %s --help' % clientName
   print


#---------------------------------------------------------------------------

def psm(inmatches,outak,bandid,niter,thresholdit,ksolve,bsolve):

   # Some values which need to be passed either in the psm argument list
   # or derived from info from the inmatches file...
   # These are part of the PSMFIT table schema, and mjdlo and mjdhi are 
   # also used for a QA plot.
   nite = '20130221'               # nite of observation
   run = '20130523114356_20130221' # processing run name
   project = 'OPS'                 # project name
   mjdlo = -1.                     # time and date of earliest observation in fit (expressed as an MJD) 
   mjdhi = -1.                     # time and date of latest observation in fit (expressed as an MJD) 
   psmfit_id_last = 93234          # current largest psmfit_id in the PSMFIT table
   psmversion = 'pyPSM_v0.1'       # psm version id
   mag_type = 'mag_psf'            # type of magnitude used in fit (e.g., mag_psf, mag_aper_8, ...)
   
   # Currently, pyPSM always solves for the a coefficients.
   # If, in the future, the option is added to choose whether
   # or not the a coefficients will be solved for, an "--asolve"
   # toggle should be added to Main method...
   asolve = 1

   # Extract info from the list of arguments passed to psm.py...
   # First, convert certain arguments from strings to numbers...
   bandid = int(bandid)
   niter = int(niter)
   thresholdit = float(thresholdit)

   # Determine which filter band is to be solved for
   bandlist = ['u','g','r','i','z','Y']
   print len(bandlist)
   if bandid < len(bandlist):
      fitband = bandlist[bandid]
   else:
      print 'Bandid %d has no corresponding filter...' % bandid
      sys.exit(1)
   #endif
   
   # Set certain defaults for each filter...
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

   # The default zeropoint...
   defaultzp = 0.0

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
      fd=open(infile)
      fd.readline()

      # Read input file line-by-line and populate AA and BB matrices 
      # (to solve matrix equation AA.XX=BB)...
      for l in fd:
         lsp=l.strip().split(',')
         band=lsp[18]
         if band != fitband:
            continue
         #endif
 
         colorstd=float(lsp[4])
         if ((band == 'i' or band == 'z') and (colorstd < 0.0 or colorstd > 0.7)) or ((band == 'g' or band == 'r') and (colorstd < 0.2 or colorstd > 1.2)):
            continue
         #endif

         ccd = int(lsp[10])
         magstd = float(lsp[3])
         mag = float(lsp[14])
         magerr = float(lsp[15])
         #totMagErr = math.sqrt(baseMagErr*baseMagErr + magerr*mag)
         totMagErr = baseMagErr
         weight = 1.0/(baseMagErr*baseMagErr)
         exptime = float(lsp[11])
         airmass = float(lsp[12])
         dm  = (mag - defaultzp + 2.5*math.log(exptime,10)) - magstd

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
         BB[iparam_a] = -9999.00+25.00
         # Since we always solve for a (asolve always equals 1), 
         # any CCD with nstar=0 reduces the number of free 
         # parameters by 1...
         if asolve == 1:
            nFreeParam += -1
         #endif

         # Fix b's...
         AA[iparam_b,:] = 0.0
         AA[iparam_b][iparam_b] = 1.0
         BB[iparam_b] = -9999.00
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

      # Output to matched catalog output file those stars that
      # survive this iteration's clipping...
      ofd3=open(resfile,'w')
      ofd3.write('res,airmass,magstd,colorstd,ccd,expnum\n')
      ofd2=open(outfile,'w')
      fd=open(infile)
      for l in fd:
         lsp=l.strip().split(',')
         band=lsp[18]
         if band != fitband:
            continue

         colorstd=float(lsp[4])
         if ((band == 'i' or band == 'z') and (colorstd < 0.0 or colorstd > 0.7)) or ((band == 'g' or band == 'r') and (colorstd < 0.2 or colorstd > 1.2)):
            continue

         ccd=int(lsp[10])
         magstd=float(lsp[3])
         mag=float(lsp[14])
         exptime=float(lsp[11])
         airmass=float(lsp[12])
         expnum=int(lsp[9])
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
            ofd2.write(l)
            resOutputLine = '%.4f,%.3f,%.4f,%.4f,%d,%d\n' % (dm,airmass,magstd,colorstd,ccd,expnum)
            ofd3.write(resOutputLine)
         #endif
   
      #endfor (fd)

      ofd3.close()
      ofd2.close()

      infile = outfile
      outfile = infile+'.tmp'

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
         if verbose > 1:
            print 'avg:'+str(avg)+'  avgchi:'+str(avgchi)+'  chisq:'+str(chisq)
            print 'k:'+str(XX[0])+' rms:'+str(sigma)+' n:'+str(ninfit)
            print 'kerr:'+str(math.sqrt(AAinv[0][0]))
         #endif
      else:
         print 'Number of free parameters ('+str(nFreeParam)+') >= number of data points ('+str(ninfit)+')'
         print 'Exiting now!'
         sys.exit(1)
      #endif

   #endfor (iiter)

   # Output the results of fit...
   os.system("/bin/rm -f "+inmatches+'.'+fitband+'.tmp*')
   ofd=open(outak,'w')
   ofd.write('niter_'+band+' '+str('%d'%int(niter))+'\n')
   ofd.write('thresholdit_'+band+' '+str('%.4f'%float(thresholdit))+'\n')
   ofd.write('n_'+band+' '+str('%d'%int(ninfit))+'\n')
   ofd.write('dof_'+band+' '+str('%d'%int(dof))+'\n')
   ofd.write('rms_'+band+' '+str('%.4f'%float(sigma))+'\n')
   ofd.write('chisq_'+band+' '+str('%.4f'%float(chisq))+'\n')
   ofd.write('k_'+band+' '+str('%.3f'%float(XX[0]))+' +/- '+str('%.3f'%float(math.sqrt(AAinv[0][0])))+'\n')
   for i in range(1,63):
      outputLine = 'a_%s %2d %.3f %.3f \n' % (band, i, XX[i]-25., math.sqrt(AAinv[i][i]))
      ofd.write('a_'+band+' '+str(i)+' '+str('%.3f'%float(XX[i]-25.))+'\n')
      if verbose > 1: print outputLine,
   #endfor
   for i in range(1,63):
      outputLine = 'a_%s %2d %.3f %.3f \n' % (band, i, XX[nccd+i], math.sqrt(AAinv[nccd+i][nccd+i]))
      ofd.write('b_'+band+' '+str(i)+' '+str('%.3f'%float(XX[nccd+i]))+'\n')
      if verbose > 1: print outputLine,
   #endfor
   ofd.close()

   # Output the results of fit into a FITS table file...
   os.system("/bin/rm -f "+outfitsfile)
   ccdidArray = numpy.arange(1,63)
   psmfit_idArray = ccdidArray+psmfit_id_last
   niteFormat = 'A%d' % len(nite)
   niteArray = numpy.array([nite]*62)
   mjdloArray = mjdlo*numpy.ones(nccd,dtype=numpy.int)
   mjdhiArray = mjdhi*numpy.ones(nccd,dtype=numpy.int)
   ccdidArray = numpy.arange(1,63)
   filterFormat = 'A%d' % len(band)
   filterArray = numpy.array([band]*62)
   # Note:  verify that these errors are what we expect:
   errors = numpy.sqrt(numpy.diagonal(AAinv))
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
   cfilterFormat = 'A%d' % len(cfitband)
   cfilterArray = numpy.array([cfitband]*62)
   stdcolor0Array = color0*numpy.ones(nccd,dtype=numpy.int)
   asolveArray = asolve*numpy.ones(nccd,dtype=numpy.int)
   bsolveArray = bsolve*numpy.ones(nccd,dtype=numpy.int)
   ksolveArray = ksolve*numpy.ones(nccd,dtype=numpy.int)
   runFormat = 'A%d' % len(run)
   runArray = numpy.array([run]*62)
   projectFormat = 'A%d' % len(project)
   projectArray = numpy.array([project]*62)
   mag_typeFormat = 'A%d' % len(mag_type)
   mag_typeArray = numpy.array([mag_type]*62)
   
   col1  = pyfits.Column(name='psmfit_id',      format='J', array=psmfit_idArray)
   col2  = pyfits.Column(name='nite',           format=niteFormat, array=niteArray)
   col3  = pyfits.Column(name='mjdlo',          format='E', array=mjdloArray)
   col4  = pyfits.Column(name='mjdhi',          format='E', array=mjdhiArray)
   col5  = pyfits.Column(name='ccdid',          format='I', array=ccdidArray)
   col6  = pyfits.Column(name='filter',         format=filterFormat, array=filterArray)
   col7  = pyfits.Column(name='a',              format='E', array=aArray)
   col8  = pyfits.Column(name='aerr',           format='E', array=aerrArray)
   col9  = pyfits.Column(name='b',              format='E', array=bArray)
   col10 = pyfits.Column(name='berr',           format='E', array=berrArray)
   col11 = pyfits.Column(name='k',              format='E', array=kArray)
   col12 = pyfits.Column(name='kerr',           format='E', array=kerrArray)
   col13 = pyfits.Column(name='rms',            format='E', array=rmsArray)
   col14 = pyfits.Column(name='chi2',           format='E', array=chisqArray)
   col15 = pyfits.Column(name='dof',            format='I', array=dofArray)
   col16 = pyfits.Column(name='photomtricflag', format='I', array=photometricFlagArray)
   col17 = pyfits.Column(name='psmversion',     format=psmversionFormat, array=psmversionArray)
   col18 = pyfits.Column(name='fit_timestamp',  format=fit_timestampFormat, array=fit_timestampArray)
   col19 = pyfits.Column(name='cfilter',        format=cfilterFormat, array=cfilterArray)
   col20 = pyfits.Column(name='stdcolor0',      format='E', array=stdcolor0Array)
   col21 = pyfits.Column(name='asolve',         format='I', array=asolveArray)
   col22 = pyfits.Column(name='bsolve',         format='I', array=bsolveArray)
   col23 = pyfits.Column(name='ksolve',         format='I', array=ksolveArray)
   col24 = pyfits.Column(name='run',            format=runFormat, array=runArray)
   col25 = pyfits.Column(name='project',        format=projectFormat, array=projectArray)
   col26 = pyfits.Column(name='mag_type',       format=mag_typeFormat, array=mag_typeArray)

   cols=pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21, col22, col23, col24, col25, col26])
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
      opts,args = getopt.getopt(sys.argv[6:],'h',['bsolve', 'ksolve', 'help', 'verbose='])
   except getopt.GetoptError:
      usage()
      sys.exit(1)
   #end try

   verbose = 0     # Default value for verbosity
   bsolve  = 0     # Default value for bsolve (0=no, 1=yes)
   ksolve  = 0     # Default value for ksolve (0=no, 1=yes)

   for o, a in opts:
      if o == '--bsolve':
         bsolve = 1
      elif o == '--ksolve':
         ksolve = 1
      elif o in ('-h', '--help'):
         usage()
         sys.exit(0)
      elif o in ('--verbose'):
         verbose = int(a)
      #endif
   #endfor

   # Call psm method
   psm(inmatches,outak,bandid,niter,thresholdit,ksolve,bsolve)


#---------------------------------------------------------------------------


