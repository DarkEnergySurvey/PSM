#!/usr/bin/env python
import math
import sys
import os
import re

def advanceheadlist(fid,radegl,decdegl,magul,maggl,magrl,magil,magzl,magyl,ra2,radeg1,fd1,fd2,done1,radct,decdct,expnumdct,ccddct,exptimedct,airmassdct,ampdct,magdct,magerrdct,class_stardct,flagsdct,banddct):
 done2 = 0
 if (done1==1):
  return (done1,radeg1)
 tol=float(3.0/3600.)
 delta=float(radeg1-ra2)
 if (delta > tol):
  return (done1,radeg1)
 while done2 == 0:
  l1=fd1.readline()
  if l1 ==  "":
   return (1,radeg1)
  else:
   l1s=l1.strip().split(',')
   radeg0=float(l1s[0])
   #print "should I add:"+' '+str(radeg0)+' ra2:'+str(ra2)+' '+str(radeg0-ra2-tol)
   if (radeg0-ra2) < -tol:
    return (done1,radeg1)

   radegl.append(radeg0)
   decdegl.append(float(l1s[1]))
   magul.append(float(l1s[2]))
   maggl.append(float(l1s[3]))
   magrl.append(float(l1s[4]))
   magil.append(float(l1s[5]))
   magzl.append(float(l1s[6]))
   magyl.append(float(l1s[7]))
   radct.append([])
   decdct.append([])
   expnumdct.append([])
   ccddct.append([])
   exptimedct.append([])
   airmassdct.append([])
   ampdct.append([])
   magdct.append([])
   magerrdct.append([])
   class_stardct.append([])
   flagsdct.append([])
   banddct.append([])
   #print "ral:"+str(len(radegl))+' '+str(radegl)+' ra2:'+str(ra2)
   #print "decl:"+str(decl)
   delta = (radeg1-radeg0)   
   if (delta > tol):
     done2 = 1
 
  return (done1,radeg0)

def cleancurlist(fid,radegl,decdegl,magul,maggl,magrl,magil,magzl,magyl,radct,decdct,expnumdct,ccddct,exptimedct,airmassdct,ampdct,magdct,magerrdct,class_stardct,flagsdct,banddct,ofd,ccdcount):
  #print "cleaning check for "+str(len(radegl))+' '+str(radegl)+' '+str(len(radct[0]))+' '+str(radct[0])
  i=0
  #print " rai "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(radct[i])
  #print " deci "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(decdct[i])
  #print " expnumi "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(expnumdct[i])
  #print " ccdi "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(ccddct[i])
  #print " exptimei "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(exptimedct[i])
  #print " airmassi "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(airmassdct[i])
  #print " ampi "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(ampdct[i])
  #print " magi "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(magdct[i])
  #print " magerri "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(magerrdct[i])
  #print " classstari "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(class_stardct[i])
  #print " flagsi "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(flagsdct[i])
  #print " bandi "+' '+str(i)+"check j "+str(len(radct[i]))+' '+str(banddct[i])
  if len(radct[i]) > 0:
   fid += 1
  for j in range(0,len(radct[i])):
   #print ral[i],decl[i],fidl[i],len(radct[i])
   band = banddct[i][j]
   stdmag = 0
   stdcolor = 0
   bandid = -1
   if band == 'u':
    bandid = 0
    stdmag=magul[i]
    stdcolor=magul[i]-maggl[i]
   if band == 'g':
    bandid = 1
    stdmag=maggl[i]
    stdcolor=maggl[i]-magrl[i]
   if band == 'r':
    bandid = 2
    stdmag=magrl[i]
    stdcolor=maggl[i]-magrl[i]
   if band == 'i':
    bandid = 3
    stdmag=magil[i]
    stdcolor=magil[i]-magzl[i]
   if band == 'z':
    bandid = 4
    stdmag=magzl[i]
    stdcolor=magil[i]-magzl[i]
   if band == 'y' or band == 'Y':
    bandid = 5
    stdmag=magyl[i]
    stdcolor=magzl[i]-magyl[i]
   ccdcount[ccddct[i][j]] += 1
   ofd.write(str(fid)+','+str(radegl[i])+','+str(decdegl[i])+','+str(stdmag)+','+str(stdcolor)+','+str(bandid)+','+str(j)+','+str('%.6f'%radct[i][j])+','+str('%.6f'%decdct[i][j])+','+str(expnumdct[i][j])+','+str(ccddct[i][j])+','+str(exptimedct[i][j])+','+str(airmassdct[i][j])+','+str(ampdct[i][j])+','+str(magdct[i][j])+','+str('%.4f'%magerrdct[i][j])+','+str(class_stardct[i][j])+','+str(flagsdct[i][j])+','+str(banddct[i][j])+'\n')

  del radct[i]
  del decdct[i]
  del expnumdct[i]
  del ccddct[i]
  del exptimedct[i]
  del airmassdct[i]
  del ampdct[i]
  del magdct[i]
  del magerrdct[i]
  del class_stardct[i]
  del flagsdct[i]
  del banddct[i]

  del radegl[i]
  del decdegl[i]
  del magul[i]
  del maggl[i]
  del magrl[i]
  del magil[i]
  del magzl[i]
  del magyl[i]
  return fid

  
def incurlist(radegl,decdegl,ra2,dec2,expnum2,ccd2,exptime2,airmass2,amp2,mag2,magerr2,class_star2,flags2,band2,radct,decdct,expnumdct,ccddct,exptimedct,airmassdct,ampdct,magdct,magerrdct,class_stardct,flagsdct,banddct):
  dtr=float(3.1415926/180.0)
  tol=float((2.0/3600.)**2)
  cosd=math.cos(dec2*dtr)
  for i in range(0,len(radegl)):
   #print 'i is:'+str(i)+' ral is'+str(radegl)
   delta=(ra2-radegl[i])*(ra2-radegl[i])*cosd*cosd+(dec2-decdegl[i])*(dec2-decdegl[i])
   #print str(ra2)+' '+str(dec2)+' '+str(delta)+' '+str(tol)+' '+str(float(delta)-float(tol))
   if float(delta) < float(tol):
    #print 'here adding'+' '+str(ra2)+' '+str(delta)+' '+str(tol)
    radct[i].append(ra2)
    decdct[i].append(dec2)
    expnumdct[i].append(expnum2)
    ccddct[i].append(ccd2)
    exptimedct[i].append(exptime2)
    airmassdct[i].append(airmass2)
    ampdct[i].append(amp2)
    magdct[i].append(mag2)
    magerrdct[i].append(magerr2)
    class_stardct[i].append(class_star2)
    flagsdct[i].append(flags2)
    banddct[i].append(band2)
    #print "appended another continuation of radeg1:"+str(radeg1)+" radct is:"+str(radct)
    break
   #else:
    #print "no match"
  
  return 0
  
def Usage():
	sys.exit("matchsort: file1 file2 outfile")

def matchsort(f1,f2,outfile):
 #print f1,f2,outfile

 ccdcount=[]
 for i in range(0,63):
  ccdcount.append(0)

 radct=[]
 decdct=[]
 expnumdct=[]
 ccddct=[]
 exptimedct=[]
 airmassdct=[]
 ampdct=[]
 magdct=[]
 magerrdct=[]
 class_stardct=[]
 flagsdct=[]
 banddct=[]

 radegl=[]
 decdegl=[]
 magul=[]
 maggl=[]
 magrl=[]
 magil=[]
 magzl=[]
 magyl=[]
 fidl=[]

 ofd=open(outfile,'w')
 ofd.write('fid,radegstd,decdegstd,magstd,colorstd,bandidstd,j,ra,dec,expnum,ccd,exptime,airmass,amp,mag,magerr,class_star,flags,band\n')
 fid=0
 fd1=open(f1)
 h1=fd1.readline()
 h1n=h1.strip().split(',')
 for i in range(0,len(h1n)):
   #print h1n[i],i
   if h1n[i] == 'radeg':
    radegcol1=i
   if h1n[i] == 'decdeg':
    decdegcol1=i
   if h1n[i] == 'magu':
    magucol1=i
   if h1n[i] == 'magg':
    maggcol1=i
   if h1n[i] == 'magr':
    magrcol1=i
   if h1n[i] == 'magi':
    magicol1=i
   if h1n[i] == 'magz':
    magzcol1=i
   if h1n[i] == 'magy':
    magycol1=i
    
 
 fd2=open(f2)
 h2=fd2.readline()
 h2n=h2.strip().split(',')
 for i in range(0,len(h2n)):
   print h2n[i],i
   if h2n[i] == 'ra':
    racol2=i
   if h2n[i] == 'dec':
    deccol2=i
   if h2n[i] == 'expnum':
    expnumcol2=i
   if h2n[i] == 'ccd':
    ccdcol2=i
   if h2n[i] == 'exptime':
    exptimecol2=i
   if h2n[i] == 'airmass':
    airmasscol2=i
   if h2n[i] == 'amp':
    ampcol2=i
   if h2n[i] == 'mag':
    magcol2=i
   if h2n[i] == 'magerr':
    magerrcol2=i
   if h2n[i] == 'class_star':
    class_starcol2=i
   if h2n[i] == 'flags':
    flagscol2=i
   if h2n[i] == 'band':
    bandcol2=i
   
 done=0
 done1=0
 radeg1=-999
 decdeg1=-999
 ra2=-999
 dec2=-999
 tol = 3/3600.0
 linecnt = 0
 while not done:
  linecnt += 1
  if linecnt/1000.0 == int(linecnt/1000.0): 
    print linecnt
  l2=fd2.readline()
  if l2 == "":
   done = 1
   break
  else:
   l2s=l2.strip().split(',')
   ra2=float(l2s[racol2])
   dec2=float(l2s[deccol2])
   expnum2=int(l2s[expnumcol2])
   ccd2=int(l2s[ccdcol2])
   exptime2=float(l2s[exptimecol2])
   airmass2=float(l2s[airmasscol2])
   amp2=int(l2s[ampcol2])
   mag2=float(l2s[magcol2])
   magerr2=100.0*float(l2s[magerrcol2])
   class_star2=float(l2s[class_starcol2])
   flags2=int(l2s[flagscol2])
   band2=str(l2s[bandcol2])

  doneheadadv=0
  while not doneheadadv:
   (done1,radeg1)=advanceheadlist(fid,radegl,decdegl,magul,maggl,magrl,magil,magzl,magyl,ra2,radeg1,fd1,fd2,done1,radct,decdct,expnumdct,ccddct,exptimedct,airmassdct,ampdct,magdct,magerrdct,class_stardct,flagsdct,banddct)
   if (done1 == 1):
    break
   if (radeg1-ra2 > tol):
    break
   #break


  incurlist(radegl,decdegl,ra2,dec2,expnum2,ccd2,exptime2,airmass2,amp2,mag2,magerr2,class_star2,flags2,band2,radct,decdct,expnumdct,ccddct,exptimedct,airmassdct,ampdct,magdct,magerrdct,class_stardct,flagsdct,banddct)


  tol = float(3/3600.0)
  done3=0
  while not done3:
   if len(radegl) > 1:
    delt = ra2-radegl[0]
    #print "delt to clean is "+str(delt)+' '+str(tol)+' ra2:'+str(ra2)+' radegl0:'+str(radegl[0])
    if delt > tol:
     #print "going to clean now"
     fid=cleancurlist(fid,radegl,decdegl,magul,maggl,magrl,magil,magzl,magyl,radct,decdct,expnumdct,ccddct,exptimedct,airmassdct,ampdct,magdct,magerrdct,class_stardct,flagsdct,banddct,ofd,ccdcount)
    else:
     break
   else:
    done3 = 1
 
 done3 = 0
 while not done3:
  if len(radegl) > 0:
    fid=cleancurlist(fid,radegl,decdegl,magul,maggl,magrl,magil,magzl,magyl,radct,decdct,expnumdct,ccddct,exptimedct,airmassdct,ampdct,magdct,magerrdct,class_stardct,flagsdct,banddct,ofd,ccdcount)
  else:
    break
   
 ofd.close()
 #for i in range(0,63):
  #print str(i)+' '+str(ccdcount[i])
 sys.exit(0)
 

if __name__ == "__main__":
     if len(sys.argv[1:]) < 3 or (sys.argv[1] == '-h'):
        Usage()
     matchsort(sys.argv[1],sys.argv[2],sys.argv[3])

