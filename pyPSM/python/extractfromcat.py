#!/usr/bin/env python

import pyfits
import os
import sys
import re

#from Header:exptime filename: exposure number, ccd number, band (r)
#xwin_image
#ywin_image
#alphawin_j2000
#deltawin_j2000
#mag_psf, or flux_psf, fluxerr_psf?, magerr_psf
#
#spread_model, class_star spreaderr_model

def Usage():
   sys.exit("extractfromcat.py dataroot pipeline attempt explist(or all) ralo rahi declo dechi classstarcut flaglimit outfile")

def extractfromcat(dataroot,pipeline,attempt,explist,ralo,rahi,declo,dechi,classstarcut,flaglimit,outfile):

 ofd=open(outfile+'.tmp','w')
 if explist == "all":
  expliststring = "all"
 else:
  expliststring = []
  d=explist.split(',')
  for i in range(0,len(d)):
   expliststring.append(str(d[i]))
 print 'expliststring ',expliststring
  
 for f in os.listdir(dataroot):
  m=re.search('([0-9]+)_r'+pipeline+'p'+str('%02d'%int(attempt)),f)
  if m:
   dir=m.group(0)
   exp=m.group(1).lstrip('0')
   expnum=int(exp)

   print dir+'is dir'
   if pipeline=='first':
     predir=dataroot+'/'+dir+'/04_catlg/data'
     pattern='D%08d'%expnum+'_([a-zY])_([0-9]+)_r01p01_firstcat.fits'
   else:
     predir=dataroot+'/'+dir+'/05_flcat/data'
     pattern='D%08d'%expnum+'_([a-zY])_([0-9]+)_r01p01_fullcat.fits'

   if not os.path.isdir(predir):
    continue
 
   ifile=exp+'.info'
   exptime = -1

   fd=open(dataroot+'/'+dir+'/'+ifile)
   for line in fd:
    m=re.search('AIRMASS\s*=\s*([0-9\.EDed\-\+]+)',line)
    if m: 
     airmass = float(m.group(1))
    m=re.search('telradeg\s*=\s*([0-9\.EDed\-\+]+)',line)
    if m: 
     telradeg = float(m.group(1))
    m=re.search('teldecdeg\s*=\s*([0-9\.EDed\-\+]+)',line)
    if m: 
     teldecdeg = float(m.group(1))
    if m: 
     exptime = m.group(1)
    m=re.search('EXPTIME\s*=\s*([0-9\.]+)',line)
    if m: 
     exptime = m.group(1)

   print str(expnum)+' '+' '+str(exptime)+' '+str(telradeg)+' '+str(teldecdeg)
   if (telradeg < float(ralo) or telradeg > float(rahi) or teldecdeg < float(declo) or teldecdeg > float(dechi)):
    print "exposure out of ra dec box"
    continue
   for file in os.listdir(predir):
    m=re.search(pattern,file)
    if m:
     band=m.group(1)
     ccd = m.group(2).lstrip('0')
     if (expliststring == 'all') or (str(expnum) in expliststring):
      #print "here at ",expnum
      t=pyfits.open(predir+'/'+file)
      n=t[2].header['NAXIS2']
      d=t[2].data
      for i in range(0,int(n)):
       xpix=d.field('xwin_image')[i]
       #ypix=d.field('ywin_image')[i]
       if float(xpix) < 1024:
         amp = 0
       else:
         amp = 1
       flags=d.field('flags')[i]
       if pipeline == 'first':
        ra=d.field('alpha_j2000')[i]
        dec=d.field('delta_j2000')[i]
        mag=d.field('mag_auto')[i]
        magerr=d.field('magerr_auto')[i]
        #flux=d.field('flux_auto')[i]
        #fluxerr=d.field('fluxerr_auto')[i]
       else:
        ra=d.field('alphawin_j2000')[i]
        dec=d.field('deltawin_j2000')[i]
        #mag=d.field('mag_psf')[i]
        #magerr=d.field('magerr_psf')[i]
        mag=d.field('mag_aper')[i][7]
        magerr=d.field('magerr_aper')[i][7]
        #flux=d.field('flux_psf')[i]
        #fluxerr=d.field('fluxerr_psf')[i]
        
       #flux_psf=pow(10,-0.4*(float(d.field('mag_aper')[i][7])-25.0))
       class_star=d.field('class_star')[i]
       #spread_model=d.field('theta_image')[i]
       #spreaderr_model=d.field('errtheta_image')[i]
       #print class_star,flags,mag
       if float(class_star) > float(classstarcut) and int(flags) <= flaglimit and float(mag) > 0.0 and float(mag) < 99.0:
        ofd.write(str(expnum)+','+str(ccd)+','+str(exptime)+','+str(airmass)+','+str(amp)+','+str(ra)+','+str(dec)+','+str(mag)+','+str(magerr)+','+str(class_star)+','+str(flags)+','+str(band)+'\n')
 
 
 ofd.close()
 ofdf=open(outfile,'w')
 ofdf.write('expnum,ccd,exptime,airmass,amp,ra,dec,mag,magerr,class_star,flags,band\n')
 ofdf.close()
 os.system("sort -t, --key=6n,7n --key=2n --key=1n "+outfile+".tmp >> "+outfile)
 
 
   
 
if __name__ == "__main__":
      if len(sys.argv[1:]) < 11 or (sys.argv[1] == '-h'):
         Usage()
      extractfromcat(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10],sys.argv[11])
 
