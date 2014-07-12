#!/usr/bin/env python
"""
matchsort.py

Take a standard star catalog (pre-sorted in ascending order by RA) 
and an observed catalog (also pre-sorted in ascending order by RA) 
and match them by their RA,DEC coordinates.

Examples:

matchsort.py --help

matchsort.py --inputStdStarCatFile standard_stars_all_id6_pypsmFormat.csv --inputObsCatFile obsquery-20131002-g-r03p01.csv --outputMatchFile matched-20131002-g-r03p01.csv --verbose 1

"""

import math
import sys
import os
import re

##################################

def main():
    import argparse
    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--inputStdStarCatFile', help='CSV file containing the standard star catalog to be used in the match', default='standardstarcat.csv')
    parser.add_argument('--inputObsCatFile', help='CSV file containing the catalog of observations to be used in the match', default='obsquery.csv')
    parser.add_argument('--outputMatchFile', help='CSV file containing the output of the match', default='matched.csv')
    parser.add_argument('--verbose', help='verbosity level of output to screen (0, 1, 2, ...)', type=int, default=0)
                        
    args = parser.parse_args()

    matchsort(args)


##################################

#advanceheadlist takes an RA from the observed catalog and creates/updates
# a list of standard stars within a "sliding window" having a range of +/- 3arcsec of that RA
# and taken from a catalog of standard stars that has been pre-sorted into order of ascending RA

def advanceheadlist(fid,radegl,decdegl,magul,maggl,magrl,magil,magzl,magyl,ra2,radeg1,fd1,fd2,done1,radct,decdct,banddct,ccddct,obslinedct):

    # DLT:  I think done2 may currently be superfluous; see comments on the "while done2 == 0"
    #       block below...
    done2 = 0
    
    # if we have finished reading the observed catalog, return...
    if (done1==1):
        return (done1,radeg1)
    
    tol=float(3.0/3600.)
    
    # if the previous standard star RA (radeg1) is above the upper bound of the tolerance
    #  range, return...
    delta=float(radeg1-ra2)
    if (delta > tol):
        return (done1,radeg1)
    
    # DLT:  I don't think this "while" (or the done2 parameter) is strictly necessary,
    # since the "return (done1, radeg0)" at the end of this method occurs *within* the
    # while loop, not after it.  (Did I accidentally change the indentation of the
    # "return (done1, radeg0) in an earlier incarnation of this file???)
    while done2 == 0:
        
        # Read a line from the standard star file...
        
        l1=fd1.readline()
        
        # if we have reached the end of the standard star file,
        # return, indicating the standard star file is "done";
        # otherwise, process the new line...
        if l1 ==  "":
            return (1,radeg1)
        else:
            l1s=l1.strip().split(',')
            radeg0=float(l1s[0])
            
            # if the RA of the standard star RA (radeg0) is below the lower bound of the tolerance
            #  range, return...
            if (radeg0-ra2) < -tol:
                return (done1,radeg1)
            
            # add the standard star info to lists of ra, dec, mags for this sliding window...
            radegl.append(radeg0)
            decdegl.append(float(l1s[1]))
            magul.append(float(l1s[2]))
            maggl.append(float(l1s[3]))
            magrl.append(float(l1s[4]))
            magil.append(float(l1s[5]))
            magzl.append(float(l1s[6]))
            magyl.append(float(l1s[7]))
            
            # initialize lists for possible observed star/standard star matches and add
            # these lists to "dictionaries" associated with this standard star...
            radct.append([])
            decdct.append([])
            ccddct.append([])
            banddct.append([])
            obslinedct.append([])
            
            # if the RA of the previous standard star is above the upper bound of the
            #  tolerance range from the RA of this standard star, declare done2=1
            delta = (radeg1-radeg0)
            if (delta > tol):
                done2 = 1
        
        # DLT:  Shouldn't this return be outside the "while done2 == 0" block?
        #       With its current indentation, it is the final statement *within*
        #       the "while done2 == 0" block.
        return (done1,radeg0)

##################################

# Write out the matched standard star/observed data to matched file...

def cleancurlist(fid,radegl,decdegl,magul,maggl,magrl,magil,magzl,magyl,radct,decdct,banddct,ccddct,obslinedct,ofd,ccdcount):
    
    # we only care about entry i=0 in the "dictionaries"...  (why?)
    i=0
    
    ## if there is at least one entry in the RA dictionary for this star, increment the running star id
    #if len(radct[i]) > 0:
    #    fid += 1
    
    # Loop through all the observations matched with standard star i...
    # (Note that many standard stars may have zero matches...)
    for j in range(0,len(radct[i])):
        band = banddct[i][j]
        stdmag = 0
        stdcolor = 0
        bandid = -1
        if band == 'u':
            bandid = 0
            stdmag=magul[i]
            stdcmag=maggl[i]
            stdcolor=magul[i]-maggl[i]
        if band == 'g':
            bandid = 1
            stdmag=maggl[i]
            stdcmag=magrl[i]
            stdcolor=maggl[i]-magrl[i]
        if band == 'r':
            bandid = 2
            stdmag=magrl[i]
            stdcmag=maggl[i]
            stdcolor=maggl[i]-magrl[i]
        if band == 'i':
            bandid = 3
            stdmag=magil[i]
            stdcmag=magzl[i]
            stdcolor=magil[i]-magzl[i]
        if band == 'z':
            bandid = 4
            stdmag=magzl[i]
            stdcmag=magil[i]
            stdcolor=magil[i]-magzl[i]
        if band == 'y' or band == 'Y':
            bandid = 5
            stdmag=magyl[i]
            stdcmag=magzl[i]
            stdcolor=magzl[i]-magyl[i]

        # only include standard star in output if stdmag or stdcmag are both positive...
        if ( (stdmag>0.) and (stdcmag>0.) ): 

            # increment the running star id
            fid += 1

            ccdcount[ccddct[i][j]] += 1

            outputLine = """%d,%f,%f,%f,%f,%d,%d,%s\n""" % (fid,radegl[i],decdegl[i],stdmag,stdcolor,bandid,j,obslinedct[i][j])
            ofd.write(outputLine)



    # Delete the dictionaries associated with standard star i (=0 in the sliding RA window)...
    del radct[i]
    del decdct[i]
    del ccddct[i]
    del banddct[i]
    del obslinedct[i]

    
    # Delete the lists associated with standard star i (=0 in the sliding RA window)...
    del radegl[i]
    del decdegl[i]
    del magul[i]
    del maggl[i]
    del magrl[i]
    del magil[i]
    del magzl[i]
    del magyl[i]

    return fid

##################################

# Find first good match (not necessarily best matches) between an observed star and standard stars
#  and place info into "dictionaries"
def incurlist(radegl,decdegl,ra2,dec2,band2,ccd2,obsline2,radct,decdct,banddct,ccddct,obslinedct):

    dtr=float(3.1415926/180.0)
    tol=float((2.0/3600.)**2)
    cosd=math.cos(dec2*dtr)

    # Loop through all standards stars i in the sliding RA window for that oberved star...
    for i in range(0,len(radegl)):

        delta=(ra2-radegl[i])*(ra2-radegl[i])*cosd*cosd+(dec2-decdegl[i])*(dec2-decdegl[i])

        # Is the sky position of standard star i (in the sliding RA window) within a 2-arcsec radial tolerance
        # of the observed star?  If so, add the observed info to that standard star's dictionaries...
        if float(delta) < float(tol):
            radct[i].append(ra2)
            decdct[i].append(dec2)
            ccddct[i].append(ccd2)
            banddct[i].append(band2)
            obslinedct[i].append(obsline2)

            # if we found one match, we take it and break out of the loop...
            break
    #else: 
    #   print "no match"
    
    return 0


# The upper-level match method...
def matchsort(args):
    
    f1=args.inputStdStarCatFile
    f2=args.inputObsCatFile
    outfile=args.outputMatchFile
    verbose=args.verbose

    # initialize ccd counts for all ccds
    ccdcount=[]
    for i in range(0,63):
        ccdcount.append(0)
    
    # initialize "dictionaries".
    # Each element of "dictionary" is associated with a standard star.
    # Each element is a list of information from the potential matches from the observed data.
    radct=[]
    decdct=[]
    ccddct=[]
    banddct=[]
    obslinedct=[]
    
    # initialize lists of standard stars.
    # These are actually lists of standards within a given sliding window of RA.
    radegl=[]
    decdegl=[]
    magul=[]
    maggl=[]
    magrl=[]
    magil=[]
    magzl=[]
    magyl=[]
    
    # Open the output file for the standard star/observed star matches...
    ofd=open(outfile,'w')

    fid=0
    
    # Open the standard star CSV file...
    fd1=open(f1)
    
    # Identify columns in the standard star CSV file corresponding to
    # radeg, decdeg, magu, magg, magr, magi, magz, magy
    h1=fd1.readline()
    h1n=h1.strip().split(',')
    for i in range(0,len(h1n)):
        if h1n[i].upper() == 'RADEG':
            radegcol1=i
        if h1n[i].upper() == 'DECDEG':
            decdegcol1=i
        if h1n[i].upper() == 'MAGU':
            magucol1=i
        if h1n[i].upper() == 'MAGG':
            maggcol1=i
        if h1n[i].upper() == 'MAGR':
            magrcol1=i
        if h1n[i].upper() == 'MAGI':
            magicol1=i
        if h1n[i].upper() == 'MAGZ':
            magzcol1=i
        if h1n[i].upper() == 'MAGY':
            magycol1=i
    
    # Open CSV file of observed data...
    fd2=open(f2)
    
    # Identify columns in the CSV observed data file corresponding to ra,dec
    h2=fd2.readline()

    #   ... but first write header for the output CSV file...
    outputHeader = h2.strip().upper()
    outputHeader = """fid,radegstd,decdegstd,magstd,colorstd,bandidstd,j,%s\n""" % (outputHeader)
    ofd.write(outputHeader)

    h2n=h2.strip().split(',')
    for i in range(0,len(h2n)):
        if h2n[i].upper() == 'RADEG':
            racol2=i
        if h2n[i].upper() == 'DECDEG':
            deccol2=i
        if h2n[i].upper() == 'BAND':
            bandcol2=i
        if h2n[i].upper() == 'CCDNUM':
            ccdcol2=i


    # initialize some variables
    #  done = "are we done with the whole process yet?"
    #  done1 = "are we done reading the observations file yet?"
    #  radeg1, decdeg1 = "initial/previous values for standard star RA,DEC"
    #  ra2, dec2 = "initial/previous values for observed star RA,DEC"
    #  tol = "tolerance in arcsec"
    #  linecnt = "line count"
    done=0
    done1=0
    radeg1=-999
    decdeg1=-999
    ra2=-999
    dec2=-999
    tol = 3/3600.0
    linecnt = 0
    
    # Loop through file of observed data...
    while not done:
        linecnt += 1
        if ( (linecnt/1000.0 == int(linecnt/1000.0)) and (verbose > 1) ):
            print '\r'+'Progress (lines read from observed catalog):  ',linecnt,
            sys.stdout.flush()
        l2=fd2.readline()

        # Are we done reading through the file of observed data yet?
        # If so, break out of loop; otherwise, process the data line...
        if l2 == "":
            done = 1
            break
        else:
            #obsline2 holds the whole line of information for this entry for future use...
            obsline2=l2.strip()
            l2s=l2.strip().split(',')
            ra2=float(l2s[racol2])
            dec2=float(l2s[deccol2])
            band2=str(l2s[bandcol2])
            ccd2=int(l2s[ccdcol2])


        # Update the sliding RA window of possibly matching standard stars
        #  and find possible matches...
        doneheadadv=0
        while not doneheadadv:
            (done1,radeg1)=advanceheadlist(fid,radegl,decdegl,magul,maggl,magrl,magil,magzl,magyl,ra2,radeg1,fd1,fd2,done1,radct,decdct,banddct,ccddct,obslinedct)
            if (done1 == 1):
                break
            if (radeg1-ra2 > tol):
                break
        #break
        
        # For this observed star, fill the "dictionaries" with all good matches from the sliding window lists of possible matches
        # that fall within a radial tolerance of 2 arcsec...
        incurlist(radegl,decdegl,ra2,dec2,band2,ccd2,obsline2,radct,decdct,banddct,ccddct,obslinedct)
        
        
        tol = float(3/3600.0)
        done3=0
        while not done3:
            if len(radegl) > 1:
                delt = ra2-radegl[0]
                if delt > tol:
                    fid=cleancurlist(fid,radegl,decdegl,magul,maggl,magrl,magil,magzl,magyl,radct,decdct,banddct,ccddct,obslinedct,ofd,ccdcount)
                else:
                    break
            else:
                done3 = 1
    
    done3 = 0
    while not done3:
        if len(radegl) > 0:
            fid=cleancurlist(fid,radegl,decdegl,magul,maggl,magrl,magil,magzl,magyl,radct,decdct,banddct,ccddct,obslinedct,ofd,ccdcount)
        else:
            break
    
    # close the input and output files...
    fd1.close()
    fd2.close()
    ofd.close()

    if verbose > 0:
        print 
        print 'ccdnum,nstds'
        for i in range(0,63):
            print str(i)+' '+str(ccdcount[i])
        print

    sys.exit(0)



##################################

if __name__ == "__main__":
    main()

##################################
