#!/usr/bin/env python

import argparse
import sys
import curses
import csv
from pathlib import Path
import numpy as np

"""
rebin.py

Rebin energy spectrum from varispaced (translated from tof spectrum) to equispaced 
energy bins.

v0.1 22/6/2020
"""

helpmessage="""
rebin.py: rebin energy spectrum from varispaced to equispaced bins.

Usage:

1) rebin.py      Run from command line with input from running terminal.
2) rebin --help  Prints this message.
3) rebin <filename> [args] 
                 Run from command line with no further input.

   args:
   -b, --binsize:    Output file binsize in keV.
                     If omitted, default 30 keV applies.
   -s, --start:      Starting energy in keV.
                     If omitted, default 0 keV applies.
   -o, --outfile:    Output filename. 
                     If omitted, input filename with binsize appended.
   -h, --help:       Print this message.
"""

help_epilog="""
There are three modes of usage:

1) rebin.py      Run from command line with input from running terminal.
2) rebin --help  Prints this message.
3) rebin [args]  <filename> 
                 Run from command line with no further input.

Input and output files are assumed to be CSV files.
"""

def processargs():
    parser = argparse.ArgumentParser(
             description="Rebin energy spectrum from varispaced to equispaced bins.",
             formatter_class=argparse.RawDescriptionHelpFormatter,
             epilog=help_epilog)
    parser.add_argument('infile', help="Input filename")
    parser.add_argument('-b','--binsize', help="Output file binsize in MeV.",
                        default=0.030, type=float)
    parser.add_argument('-s','--start', help="Starting energy in MeV", default=0.000,type=float)
    parser.add_argument('-o','--outfile', help="Output filename")
    #parser.add_argument('--', help="", default=)
    args = parser.parse_args()
    return args

def readwrite(win):
    global infile,outfile,start,binsize
    leave=b'n'
    while leave != b"y":
        win.addstr(3,0,"Enter input filename : ")
        infile=get_valid_name(win, default=infile)
        win.addstr(5,0,"Enter output filename: ")
        outfile=get_valid_name(win, default=outfile)
        win.addstr(7,0,"Enter binsize in MeV : ")
        binsize=get_valid_float(win, default=binsize)
        win.addstr(9,0,"Enter start energy in MeV : ")
        start=get_valid_float(win, default=start)
        win.addstr(15,0,"OK to proceed? [y/n]: ")
        y, x =win.getyx()
        curses.echo()
        leave=win.getstr(y,x)
        curses.noecho()
        #S2=str(s2,encoding='utf-8')

def get_valid_name(win, default=None):
    y, x =win.getyx()
    curses.echo()
    s=str(win.getstr(y,x,30),encoding='utf-8')
    if s=="": s=default
    curses.noecho()
    return s

def get_valid_float(win, default=None):
    y, x =win.getyx()
    curses.echo()
    loop=True
    while loop:
        try:
            if default != None: win.addstr(y,x,"%.3f"%default)
            s=str(win.getstr(y,x,10),encoding='utf-8')
            if s=="":
                f=default
            else:
                f=float(s)
            loop=False
        except:
            pass
    curses.noecho()
    return f
    
class boundary(object):
    """
    Utility class to hold details of histogram bin
    """
    def __init__(self, index, count, low, high):
        self.index=index  # index in input array
        self.count=count  # value of count (integer) in bin
        self.low=low      # low edge of bin
        self.high=high    # upper edge of bin

def GetBorders(abscissa, counts):
    """
    Gives the lower and upper bound of a channel in the spectrum
    input:   e   numpy array containing values for abscissa
    """
    if not isinstance(abscissa, np.ndarray): raise ValueError("Input 'abscissa' is not an array")
    N=len(abscissa)
    if N<2: raise ValueError("Too few points in abscissa array")
    if not isinstance(counts, np.ndarray): raise ValueError("Input 'counts' is not an array")
    Nc=len(counts)
    if Nc != N: raise ValueError("Input arrays have different lengths")
            
    bounds=[]
        
    for i in range(N):
        index=i
        if i==0: 
            dup=(abscissa[i+1]-abscissa[i])/2
            low = abscissa[i]-dup
            high = abscissa[i]+dup
            b = boundary(index, counts[i], low, high)
        elif i==N-1:
            ddn=(abscissa[i]-abscissa[i-1])/2
            low = abscissa[i]-ddn
            high = abscissa[i]+ddn
            b = boundary(index, counts[i], low, high)
        else:
            dup=(abscissa[i+1]-abscissa[i])/2
            ddn=(abscissa[i]-abscissa[i-1])/2
            low = abscissa[i]-ddn
            high = abscissa[i]+dup
            b = boundary(index, counts[i], low, high)
        bounds.append(b)
            
    return bounds

def SplitCount(count,inspan,outspan):
        """
        Split count into channel count + leftover falling outside channel
        """
        if count == 0: return (0, 0)
        inwidth = inspan[1]-inspan[0]
        outwidth = outspan[1]-outspan[0]
        if inwidth < outwidth: return (count,0)
        fraction=outwidth/inwidth
        if fraction <0.0 or fraction >1.0:
            print("Fraction error")
        fuzz=0.5
        newcount = int(fraction*count+fuzz)
        leftover = count-newcount
        return newcount, leftover

def Rebinner(old, counts, new, debug=False):
    """
    Rebins spectrum
    input:   old     numpy array containing values for original scale
             count   numpy array containing values for original counts
             new     numpy array containing values for new scale, equally spaced
    """
    DEBUG=debug
    if not isinstance(old, np.ndarray): raise ValueError("Input old is not an array")
    if not isinstance(new, np.ndarray): raise ValueError("Input new is not an array")
    if not isinstance(counts, np.ndarray): raise ValueError("Input counts is not an array")
    N1=len(old)
    if N1<2: raise ValueError("Too few points in original scale array")
    N2=len(counts)
    print(N1,N2)
    if N2<2: raise ValueError("Too few points in original count array")
    if N1 != N2: raise ValueError("old and count do not have same length")
    N3=len(new)
    if N3<2: raise ValueError("Too few points in new scale array")
        
    newcount=np.zeros(N3)

    inbounds=GetBorders(old, counts)
    outdiff=(new[1]-new[0])/2
    inlowest=inbounds[0].low
    inhighest=inbounds[-1].high
    inew=0
    iold=0
    leftover=0
    
    # skip if below lower limit in input
    while new[inew]+outdiff < inlowest:
        newcount[inew]=0
        inew+=1
        
    # fill output bin by bin
    for i in range(inew,N3):
        # borders of output bin
        nlow=new[i]-outdiff
        nhigh=new[i]+outdiff
        if DEBUG: print('infor',i,'target:',nlow,nhigh,'source',inbounds[iold].low,inbounds[iold].high)
        # skip if above upper limit in input
        if nlow > inhighest:
            newcount[i]=0
            continue
        # if input bins are narrow keep incrementing same output bin until border reached
        while iold<N1-1 and inbounds[iold].high < nhigh:
            newcount[i] += inbounds[iold].count
            if DEBUG: print('>newcount',i,newcount[i])
            iold += 1
        # right border of input lies above right border of output: split counts 
        ospan = (inbounds[iold].low, new[i]+outdiff)
        ispan = (inbounds[iold].low, inbounds[iold].high)
        newc,leftover = SplitCount(inbounds[iold].count,ispan,ospan)
        if DEBUG: print('ispan:', ispan, ' ospan',ospan, " nc,lc:",newc,leftover)
        # increment output
        newcount[i] += newc #leftover + nc
        if DEBUG: print('newcount',i,newcount[i],leftover)
        # if input overlaps output, keep leftover counts for next output bin ...
        lo_lo=nhigh
        lo_hi=inbounds[iold].high
        if lo_hi > lo_lo:
            if DEBUG: print("skip")
            inbounds[iold].low=nhigh
            inbounds[iold].count = leftover
            leftover = 0
        else:  # ... handle next input bin
            iold+=1
        # check for out of bounds
        if iold>=N1: break
        if DEBUG: print("")
        
    return newcount

def ReadFile(name):
    with open(name,newline='') as csvfile:
        r=csv.reader(csvfile)
        data=[]
        tof=np.zeros(3100)
        count=np.zeros(3100,dtype=int)
        otchan=0
        counts=0
        for row in r:
            if len(row)<2: break
            if otchan==0: 
                otchan=int(float(row[0])*1000+0.01)
                ###print("Initial channel=", otchan)
                firstchan=otchan
            tmptof=float(row[0])
            tmpcount=int(row[1])
            tchan=int(tmptof*1000+0.01)
            tof[tchan]=tmptof
            count[tchan]+=tmpcount
    ###print("Last channel=", tchan)
    lastchan=tchan
    for i in range(3100):
        if tof[i]==0.0: continue
        data.append([tof[i],count[i]])
    data=np.array(data)
    ###print(np.shape(data))
    return data

def WriteFile(name, outE, outdata):
    if name == 'r00000': return
    with open(name,'w',newline='') as csvwfile:
        w=csv.writer(csvwfile)
        for i in range(len(outdata)):
            E="%5.3f"%(outE[i]+0.00001,)
            C="%6.0f"%(outdata[i],)
            w.writerow([E,C])

if __name__=="__main__":
 
    # globals
    start=0.0
    binsize=0.030
    infile=None
    outfile=None

    arglist=sys.argv
    print('')
    if len(arglist)>1:
        args=processargs()
        binsize=args.binsize
        start=args.start
        infile=args.infile
        outfile=args.outfile
        if outfile == None:
            outfile='r00000'
    else:
        curses.wrapper(readwrite)
    ###print("Here we go...")
    ###print(infile,outfile,start,binsize)
    Infile=Path(infile)
    if Infile.is_file(): print(Infile, "exists")
    stem=Infile.stem
    suffix=Infile.suffix
    if outfile=="" or outfile==None:
        Outfile=Path(stem+"_rebin"+suffix)
    else:
        Outfile=Path(outfile)
    if Outfile.is_file(): print(Outfile, "exists")
    ###print(Outfile)

    print("infile:",Infile)
    print("outfile:",Outfile)
    print("start:",start)
    print("binsize:",binsize)
    
    data = ReadFile(Infile)
    maxE = data[-1,0]+binsize
    print(maxE)
    outE=np.arange(start,maxE,binsize)
    print(len(outE))
    R=Rebinner(data[:,0],data[:,1],outE)
    print("R", len(R))
    WriteFile(Outfile, outE, R)
    
    #import matplotlib.pyplot as plt
    #plt.figure(figsize=(16,6))
    #scale=np.max(data[:,1])/np.max(R)
    #plt.plot(data[:,0],data[:,1]/scale,drawstyle='steps-mid')
    #plt.plot(outE,R,drawstyle='steps-mid')
    #plt.xlabel("Energy [MeV]")
    #plt.ylabel("Counts per channel")
    #plt.show()

    
    
    
