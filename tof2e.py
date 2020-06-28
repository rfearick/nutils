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

help_epilog="""
There are three modes of usage:

1) rebin.py      Run from command line with input from running terminal.
2) rebin --help  Prints this message.
3) rebin [args]  <filename> 
                 Run from command line with no further input.

Input and output files are assumed to be CSV files.
"""

def processargs():
    
    """
    Process command line arguments.
    This is a basic set of commands for Python standard lib 'argparse'
    """
    global binsize, start
    parser = argparse.ArgumentParser(
             description="Rebin energy spectrum from varispaced to equispaced bins.",
             formatter_class=argparse.RawDescriptionHelpFormatter,
             epilog=help_epilog)
    parser.add_argument('infile', help="Input filename")
    parser.add_argument('-t','--tzero', help="Time zero in channels.",
                        default=binsize, type=float)
    parser.add_argument('-d','--distance', help="Target to detector distance in m.",
                        default=binsize, type=float)
    parser.add_argument('-b','--binsize', help="Output file binsize in MeV.",
                        default=binsize, type=float)
    parser.add_argument('-c','--chperns', help="TOF channels per ns.",
                        default=binsize, type=float)
    parser.add_argument('-s','--start', help="Starting energy in MeV",
                        default=start,type=float)
    parser.add_argument('-o','--outfile', help="Output filename")
    #parser.add_argument('--', help="", default=)
    args = parser.parse_args()
    return args

def readwrite(win,t0,distance,chperns,binsize,start,Infile,Outfile):

    """
    Read rebin data from terminal.
    This uses the Python standard library 'curses'
    """
    
    leave=b'n'
    while leave != b'y':
        win.addstr(3,0,f"Enter input filename [{Infile}]: ")
        Infile=get_valid_name(win, default=Infile, mustexist=True)
        if Outfile==None: Outfile=make_new_filepath(Infile,Outfile)
        win.addstr(5,0,f"Enter output filename [{Outfile}]: ")
        Outfile=get_valid_name(win, default=Outfile)
        win.addstr(7,0,f"Enter t0 in channels [{t0:.3f}]: ")
        t0=get_valid_float(win, default=t0)
        win.addstr(9,0,f"Enter target-detector distance in m [{distance:.3f}]: ")
        distance=get_valid_float(win, default=distance)
        win.addstr(11,0,f"Enter TOF channels per ns [{chperns:.3f}]: ")
        chperns=get_valid_float(win,default=chperns)
        win.addstr(13,0,f"Enter binsize in MeV [{binsize:.3f}]: ")
        binsize=get_valid_float(win, default=binsize)
        win.addstr(15,0,f"Enter start energy in MeV [{start:.3f}]: ")
        start=get_valid_float(win, default=start)
        win.addstr(21,0,"OK to proceed? [y/n/q]: ")
        y, x =win.getyx()
        curses.echo()
        leave=win.getstr(y,x)
        curses.noecho()
        if leave==b'q': exit()
    return t0,distance,chperns,binsize, start, Infile, Outfile

def get_valid_name(win, default=None, mustexist=False):
    """
    Helper function for readwrite().
    Validate string.
    """
    y, x =win.getyx()
    curses.echo()
    while 1:
        win.clrtoeol()
        s=str(win.getstr(y,x,80),encoding='utf-8')
        if s=="": s=str(default)
        p=Path(s)
        if not mustexist: break
        if p.is_file(): break
        else:
            win.addstr(y+1,0,"File does not exist: enter valid name")
            win.move(y,x)
    win.move(y+1,0)
    win.clrtoeol()        
    curses.noecho()
    return p

def make_new_filepath(infile,outfile):
    """
    Helper function for readwrite.
    Create new file path for Outfile.
    """
    if outfile=="" or outfile==None:
        parent=infile.parent
        stem=infile.stem
        suffix=infile.suffix
        newname=stem+"_rebin"+suffix
        if parent==".":
            Outfile=Path(newname)
        else:
            Outfile=Path(parent)/newname
    else:
        Outfile=Path(outfile)
    return Outfile

def get_valid_float(win, default=None):
    """
    Helper function for readwrite().
    Validate float.
    """
    y, x =win.getyx()
    win.clrtoeol()
    curses.echo()
    loop=True
    while loop:
        try:
            s=str(win.getstr(y,x,20),encoding='utf-8')
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
    input:   abscissa   numpy array containing values for abscissa
             counts     numpy array containing spectrum counts
    """
    
    if not isinstance(abscissa, np.ndarray):
        raise ValueError("Input 'abscissa' is not an array")
    N=len(abscissa)
    if N<2:
        raise ValueError("Too few points in abscissa array")
    if not isinstance(counts, np.ndarray):
        raise ValueError("Input 'counts' is not an array")
    Nc=len(counts)
    if Nc != N:
        raise ValueError("Input arrays have different lengths")
            
    bounds=[]
        
    for i in range(N):
        index=i
        # Put boundaries midway between spectrum channels.
        # For first and last take both boundaries from the one that can be computed.
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
        Split count into channel count + leftover falling outside channel.
        
        input:   count        value to be split into new bin count + leftover
                 inspan       tuple of lower and upper limits of bin in old data
                 outspan      tuple of lower and upper limits of bin in new data

        The count is an integer and count = newcount+leftover holds.
        This version just rounds the integers, but fuzz can be replaced with a
        uniform random number for a perhaps better visual effect.
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
    Rebins spectrum.
    input:   old     numpy array containing values for original scale
             count   numpy array containing values for original counts
             new     numpy array containing values for new scale, equally spaced
    """
    DEBUG=debug
    if not isinstance(old, np.ndarray):
        raise ValueError("Input old is not an array")
    if not isinstance(new, np.ndarray):
        raise ValueError("Input new is not an array")
    if not isinstance(counts, np.ndarray):
        raise ValueError("Input counts is not an array")
    N1=len(old)
    if N1<2:
        raise ValueError("Too few points in original scale array")
    N2=len(counts)
    if N2<2:
        raise ValueError("Too few points in original count array")
    if N1 != N2:
        raise ValueError("old and count do not have same length")
    N3=len(new)
    if N3<2:
        raise ValueError("Too few points in new scale array")
        
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
        if DEBUG: print('infor',i,'target:',nlow,nhigh,'source',
                        inbounds[iold].low,inbounds[iold].high)
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
    """
    Read csv file containing one columns of data (count).

    Conflate all bins with same energy (due to rounding issues at source)
    while preserving total count in spectrum.
    Returns a numpy array.
    """
    with open(name,newline='') as csvfile:
        r=csv.reader(csvfile)
        data=[]
        for row in r:
            if len(row)<1: break
            #print(row)
            tmpt=float(row[0])
            tmpcount=int(row[1])
            data.append([tmpt,tmpcount])
    data=np.array(data)
    return data

def WriteFile(name, outE, outdata):
    """
    Write a csv file containing the new data in two columns.
    """
    #if name == 'r00000': return
    with open(name,'w',newline='') as csvwfile:
        w=csv.writer(csvwfile)
        for i in range(len(outdata)):
            E="%5.3f"%(outE[i]+0.00001,)
            C="%6.0f"%(outdata[i],)
            w.writerow([E,C])

if __name__=="__main__":

    c=0.299792458 # m/ns
    mn=939.565 # MeV/c^2
    # globals
    t0=0.0       # time of flight zero (time at target)
    distance=1.0 # targer to detector distance
    chperns=1.0  # channels per ns
    start=0.0    # initial energy in energy spectrum
    binsize=0.030 # width of bin in energy spectrum
    #infile=None
    #outfile=None
    Infile=None
    Outfile=None

    arglist=sys.argv
    print('')
    if len(arglist)>1:
        args=processargs()
        binsize=args.binsize
        start=args.start
        infile=args.infile
        outfile=args.outfile
        Infile=Path(infile)
        if Infile.is_file(): print(Infile, "exists")
        if outfile=="" or outfile==None:
            Outfile=make_new_filepath(Infile,outfile)
        else:
            Outfile=Path(outfile)
        if Outfile.is_file(): print(Outfile, "exists")
    else:
        ret=curses.wrapper(readwrite,t0,distance,chperns,binsize,start,Infile,Outfile)
        t0,distance,chperns,binsize,start,Infile,Outfile=ret

    print("infile:",Infile)
    print("outfile:",Outfile)
    print("t0:",t0)
    print("distance:",distance)
    print("chperns:",chperns)
    print("start:",start)
    print("binsize:",binsize)
    
    data = ReadFile(Infile)
    # extract columns of count
    tn=np.arange(len(data[:,1]))/chperns
    cn=data[:,1]
    nt0=int(t0-chperns*(distance/c))-2
    print("nt0:",nt0)
    eb=66.0*1.15
    gammab=eb/mn+1.0
    taub=gammab/np.sqrt(gammab**2-1.0)
    ntb=int(taub*(distance/c)*chperns)
    print("ntb",ntb)
    nt0=int(t0-ntb)
    cn=cn[0:nt0]
    tn=tn[0:nt0]
    #cn=cn[0:720]
    #tn=tn[0:720]
    #print(tn)
    # dimensionless time
    t0=t0/chperns
    taun=(t0-tn)/(distance/c)
    n=len(taun)-1
    while n>1:
        if taun[n]>1.0: break
        n=n-1
    taun=taun[0:n]
    cn=cn[0:n]
    #print(taun)
    #print(cn)
    # convert to gamma
    gamman=taun/np.sqrt(taun**2-1)
    # convert to kinetic energy
    En=mn*(gamman-1)
    #print(En)
    # set up target energies
    maxE = En[-1]+binsize
    print(start,maxE,binsize,En[-1])
    outE=np.arange(start,maxE,binsize)
    # rebin and write out.
    print(len(tn),len(cn),len(taun),len(gamman),len(En),len(outE))
    R=Rebinner(En,cn,outE)
    WriteFile(Outfile, outE, R)

    Plot=False
    Plot=True
    if Plot:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(16,6))
        scale=np.max(cn)/np.max(R)
        plt.plot(outE,R,drawstyle='steps-mid',label='rebinned')
        plt.plot(En,cn/scale,drawstyle='steps-mid',label='from tof')
        plt.xlabel("Energy [MeV]")
        plt.ylabel("Counts per channel")
        plt.legend()
        plt.show()

    
    
    
