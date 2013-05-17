#!/usr/bin/env python

import pyfits
import os
import sys
import re


def Usage():
   sys.exit("extractfromstdcat.py catalogfile ralo rahi declo dechi outfile")

def extractfromstdcat(catalogfile,ralo,rahi,declo,dechi,outfile):

 ofd=open(outfile+'.tmp','w')
 t=pyfits.open(catalogfile)
 n=t[1].header['NAXIS2']
 d=t[1].data
 for i in range(0,int(n)):
  radeg=float(d.field('radeg')[i])
  decdeg=float(d.field('decdeg')[i])
  if (radeg < float(ralo) or radeg > float(rahi) or decdeg < float(declo) or decdeg > float(dechi)):
   continue

  magu=d.field('stdmag_u')[i]
  magg=d.field('stdmag_g')[i]
  magr=d.field('stdmag_r')[i]
  magi=d.field('stdmag_i')[i]
  magz=d.field('stdmag_z')[i]
  magy=d.field('stdmag_y')[i]
  ofd.write(str(radeg)+','+str(decdeg)+','+str(magu)+','+str(magg)+','+str(magr)+','+str(magi)+','+str(magz)+','+str(magy)+'\n')
 
 ofd.close()
 ofdf=open(outfile,'w')
 ofdf.write('radeg,decdeg,magu,magg,magr,magi,magz,magy\n')
 ofdf.close()
 os.system("sort -t, --key=1n,2n "+outfile+".tmp >> "+outfile)
 
if __name__ == "__main__":
      if len(sys.argv[1:]) < 6 or (sys.argv[1] == '-h'):
         Usage()
      extractfromstdcat(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
 
