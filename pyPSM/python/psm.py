#!/usr/bin/env python

# Authors:   Brian Yanny and Douglas Tucker
# Date:      17 May 2013
# Updated:   18 May 2013

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

   # Solve for k coefficient?  Default is no.
   if ksolve in ['True', 'true', 'T', 't', '1', 'yes', 'y', '1']:
      ksolve = 'True'
   else:
      ksolve = 'False'
   #endif

   # Solve for b coefficients?  Default is no.
   if bsolve in ['True', 'true', 'T', 't', '1', 'yes', 'y', '1']:
      bsolve = 'True'
   else:
      bsolve = 'False'
   #endif

   # The input matched file is inmatches; 
   # the output matched file (used for iterative clipping) is inmatches.tmp; 
   # the residuals file (used for QA plots) is inmatches.res.
   infile = inmatches
   outfile = inmatches+'.'+fitband+'.tmp'
   resfile = inmatches+'.'+fitband+'.res.csv'

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
   if ksolve=='True':
      nFreeParam = nFreeParam + 1
   if bsolve=='True':
      nFreeParam = nFreeParam + nccd
   
   print ksolve, bsolve, nparam, nFreeParam

   # Iteration loop...
   for iiter in range(0,niter):
      # Initialize various variables used in the fit,
      # including the matrices AA and BB (for the
      # matrix equation AA.XX=BB)...
      sum = 0.0
      sumsq = 0.0
      ninfit = 0.0
      zv = numpy.zeros(nparam)
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
         if band == 'u':
            bdefault = 0.0
            color0 = 1.39
            kdefault = 0.489
         elif band == 'g':
            bdefault = -0.1
            color0 = 0.53
            kdefault = 0.181
         elif band == 'r':
            bdefault = -0.08
            color0 = 0.53
            kdefault = 0.095
         elif band == 'i':
            bdefault = -0.3
            color0 = 0.09
            kdefault = 0.089
         elif band == 'z':
            bdefault = -0.09
            color0 = 0.053
            kdefault = 0.089
         elif band == 'y' or band == 'Y':
            bdefault = 0.26
            color0 = 0.05
            kdefault = 0.050
         else:
            print 'Band %s is an unknown filter...' % band
            sys.exit(1)
         #endif

         ccd = int(lsp[10])
         magstd = float(lsp[3])
         mag = float(lsp[14])
         magerr = float(lsp[15])
         #weight = 1.0/(magerr+0.00001)/(magerr+0.00001)
         #weight = 1.0
         weight = 1.0/(baseMagErr*baseMagErr)
         #print mag,magerr,weight
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

            if bsolve=='True':
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

      # Set the a's, b's of any CCDs for which nstar=0 to signal values
      for iccd in range(0,len(badinds)):

         ccd = badinds[iccd] + 1
         iparam_a = ccd
         iparam_b = nccd + ccd

         # Fix a's...
         AA[iparam_a,:] = 0.0
         AA[iparam_a][iparam_a] = 1.0
         BB[iparam_a] = -9999.00+25.00

         # Fix b's...
         AA[iparam_b,:] = 0.0
         AA[iparam_b][iparam_b] = 1.0
         BB[iparam_b] = -9999.00         
      #endfor

      if (ksolve=='False'):
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
         if band == 'u':
            bdefault = 0.0
            color0 = 1.39
         if band == 'g':
            bdefault = -0.1
            color0 = 0.53
         if band == 'r':
            bdefault = -0.08
            color0 = 0.53
         if band == 'i':
            bdefault = -0.3
            color0 = 0.09
         if band == 'z':
            bdefault = -0.09
            color0 = 0.09
         if band == 'y' or band == 'Y':
            bdefault = 0.26
            color0 = 0.05
         ccd=int(lsp[10])
         magstd=float(lsp[3])
         mag=float(lsp[14])
         exptime=float(lsp[11])
         airmass=float(lsp[12])
         expnum=int(lsp[9])
         # Note:  change befault to XX[nccd+ccd]
         dm = (mag + 2.5*math.log(exptime,10) - XX[ccd] - XX[0]*airmass - bdefault*(colorstd-color0)) - magstd
         sum += dm
         sumsq += dm*dm
         ninfit += 1
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
      avg = float(sum/ninfit)
      sigma = math.sqrt(sumsq/ninfit-avg*avg)
      print "k:"+str(XX[0])+' rms:'+str(sigma)+' n:'+str(ninfit)
      print 'kerr:'+str(math.sqrt(AAinv[0][0]))

   #endfor (iiter)

   # Output the results of fit...

   os.system("/bin/rm -f "+inmatches+'.'+fitband+'.tmp*')
   ofd=open(outak,'w')
   ofd.write('n_'+band+' '+str('%d'%int(ninfit))+'\n')
   ofd.write('niter_'+band+' '+str('%d'%int(niter))+'\n')
   ofd.write('rms_'+band+' '+str('%.4f'%float(sigma))+'\n')
   ofd.write('k_'+band+' '+str('%.3f'%float(XX[0]))+'\n')
   for i in range(1,63):
      outputLine = 'a_%s %2d %.3f %.3f \n' % (band, i, XX[i]-25., math.sqrt(AAinv[i][i]))
      print outputLine,
      ofd.write('a_'+band+' '+str(i)+' '+str('%.3f'%float(XX[i]-25.))+'\n')
   for i in range(1,63):
      outputLine = 'a_%s %2d %.3f %.3f \n' % (band, i, XX[nccd+i], math.sqrt(AAinv[nccd+i][nccd+i]))
      print outputLine,
      ofd.write('b_'+band+' '+str(i)+' '+str('%.3f'%float(XX[nccd+i]))+'\n')
   ofd.close()

   print "That's all, folks!"
 
   sys.exit
 
#---------------------------------------------------------------------------
# Main method:

if __name__ == "__main__":

   # If help requested, or if insufficient number of arguments, print out usage.
   if (sys.argv[1] == '-h') or (sys.argv[1] == '--help'):
      usage()
      sys.exit(0)
   elif len(sys.argv[1:]) < 5:
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

   verbose = 0        # Default value for verbosity
   bsolve  = 'False'  # Default value for bsolve
   ksolve  = 'False'  # Default value for ksolve

   for o, a in opts:
      if o == '--bsolve':
         bsolve = 'True'
      elif o == '--ksolve':
         ksolve = 'True'
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


