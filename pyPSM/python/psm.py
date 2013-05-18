#!/usr/bin/env python

# Authors:   Brian Yanny and Douglas Tucker
# Date:      17 May 2013
# Updated:   18 May 2013

import numpy
import sys
import math
import os

#---------------------------------------------------------------------------

def Usage():
   sys.exit("psm <inmatches> <outak> <bandidugrizy012345> <niter> <thresholdit> <ksolve> <bsolve>")

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
      print bandid
      fitband = bandlist[bandid]
      print bandid, fitband
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
   # the output matched file (used for iterative clipping) is inmatches.tmp
   infile = inmatches
   outfile = inmatches+'.tmp'

   # The default zeropoint...
   defaultzp = 0.0

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
      sum = 0.0
      sumsq = 0.0
      ninfit = 0.0
      zv = numpy.zeros(nparam)
      I = numpy.array(numpy.identity(nparam))
      for i in range(0,nparam):
         for j in range(0,nparam):
            I[i][j] = 0
      b = numpy.zeros(nparam)
      #A = numpy.vstack([zv,I]).T
      A = I
      fd=open(infile)
      fd.readline()

      # Read input file line-by-line and populate A matrix and b vector (Ax=b)...
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
            bcoeff = 0.0
            color0 = 1.39
            kdefault = 0.489
         elif band == 'g':
            bcoeff = -0.1
            color0 = 0.53
            kdefault = 0.181
         elif band == 'r':
            bcoeff = -0.08
            color0 = 0.53
            kdefault = 0.095
         elif band == 'i':
            bcoeff = -0.3
            color0 = 0.09
            kdefault = 0.089
         elif band == 'z':
            bcoeff = -0.09
            color0 = 0.053
            kdefault = 0.089
         elif band == 'y' or band == 'Y':
            bcoeff = 0.26
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
         weight = 1.0
         #print mag,magerr,weight
         exptime = float(lsp[11])
         airmass = float(lsp[12])
         dm  = (mag - defaultzp + 2.5*math.log(exptime,10)) - magstd

         if abs(dm) < 4.5:

            #Code to indices:
            # kcoeff:       index 0
            # acoeff[ccd]:  index iparam_a
            # bcoeff[ccd]:  index iparam_b

            iparam_a = ccd
            iparam_b = nccd + ccd

            #Slight difference in weighting between original pyPSM version and javaPSM version.
            #Need to check which is correct.
            #Original pyPSM version:
            #b[0] += weight*weight*airmass*dm
            #b[iparam_a] += weight*dm
            #A[0][0] += weight*weight*airmass*airmass
            #A[iparam_a][0] += weight*airmass
            #A[0][iparam_a] += weight*airmass
            #A[iparam_a][iparam_a] += 1.0*weight

            #Version from javaPSM (note b[0] and A[0][0] have each lost one factor of weight):
            b[0] += weight*airmass*dm
            b[iparam_a] += weight*dm
            A[0][0] += weight*airmass*airmass
            A[iparam_a][0] += weight*airmass
            A[0][iparam_a] += weight*airmass
            A[iparam_a][iparam_a] += 1.0*weight
     
            dc = colorstd-color0

            if bsolve=='True':
               A[iparam_b][iparam_b] += dc * dc * weight
               A[0][iparam_b] += dc * airmass * weight
               A[iparam_b][0] += dc * airmass * weight
               A[iparam_a][iparam_b] += dc * weight
               A[iparam_b][iparam_a] += dc * weight
               b[iparam_b] += dc * dm * weight
            else:
               A[iparam_b][iparam_b] = 1.0
               A[0][iparam_b] += dc * airmass * weight
               A[iparam_b][0] = 0.0
               A[iparam_a][iparam_b] += dc * weight
               A[iparam_b][iparam_a] = 0.0
               #b[iparam_b] = bdefaultValues[iccd]
               b[iparam_b] = bcoeff
            #endif

         #endif

      #endfor

      if (ksolve=='False'):
         A[0][0] = 1.0
         b[0] = kdefault
         for iccd in range(0,int(nccd)):
            ccd = iccd+1
            iparam_a = ccd
            iparam_b = nccd + ccd
            A[0][iparam_a] = 0.0
            A[0][iparam_b] = 0.0
         #endfor
      #endif

      # Solve for x in Ax=b matrix equation...
      (x,residue,rank,s) = numpy.linalg.lstsq(A,b)
      print 'rank', rank

      # Output to matched catalog output file those stars that
      # survive this iteration's clipping...
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
            bcoeff = 0.0
            color0 = 1.39
         if band == 'g':
            bcoeff = -0.1
            color0 = 0.53
         if band == 'r':
            bcoeff = -0.08
            color0 = 0.53
         if band == 'i':
            bcoeff = -0.3
            color0 = 0.09
         if band == 'z':
            bcoeff = -0.09
            color0 = 0.09
         if band == 'y' or band == 'Y':
            bcoeff = 0.26
            color0 = 0.05
         ccd=int(lsp[10])
         magstd=float(lsp[3])
         mag=float(lsp[14])
         exptime=float(lsp[11])
         airmass=float(lsp[12])
         dm = (mag + 2.5*math.log(exptime,10) - x[ccd] - x[0]*airmass - bcoeff*(colorstd-color0)) - magstd
         sum += dm
         sumsq += dm*dm
         ninfit += 1
         if abs(dm) < thresholdit:
            ofd2.write(l)
         #endif
   
      #endfor (fd)

      ofd2.close()
      infile = outfile
      outfile = infile+'.tmp'
      avg = float(sum/ninfit)
      sigma = math.sqrt(sumsq/ninfit-avg*avg)
      print "k:"+str(x[0])+' rms:'+str(sigma)+' n:'+str(ninfit)
   
   #endfor (iiter)

   # Output the results of fit...

   os.system("/bin/rm -f "+inmatches+'.tmp*')
   ofd=open(outak,'w')
   ofd.write('n_'+band+' '+str('%d'%int(ninfit))+'\n')
   ofd.write('niter_'+band+' '+str('%d'%int(niter))+'\n')
   ofd.write('rms_'+band+' '+str('%.4f'%float(sigma))+'\n')
   ofd.write('k_'+band+' '+str('%.3f'%float(x[0]))+'\n')
   for i in range(1,63):
      ofd.write('a_'+band+' '+str(i)+' '+str('%.3f'%float(x[i]-25.))+'\n')
   for i in range(1,63):
      ofd.write('b_'+band+' '+str(i)+' '+str('%.3f'%float(x[nccd+i]))+'\n')
   ofd.close()

   print "That's all, folks!"
 
   sys.exit
 
#---------------------------------------------------------------------------

if __name__ == "__main__":
   if len(sys.argv[1:]) < 7 or (sys.argv[1] == '-h'):
      Usage()
   psm(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])

#---------------------------------------------------------------------------


