#from pylab import *
# leech library stuff
from ctypes import *
import os
from random import gauss,random,randint


leechdec = cdll.LoadLibrary(os.getcwd()+'/leechDecoder.so')

leechdec.decode.argtypes= [POINTER(c_float),POINTER(c_float)]
leechdec.decode.restype = c_int

floatArray =c_float*24


def dec(pts):
    cx=floatArray(*pts)
    d = c_float()
    result = leechdec.decode(cx,byref(d))

    return int(result)#,d.value
#end leech library stuff

    
def decodeD8(r):
    rInt = [0,0,0,0,0,0,0,0]
    rDist = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    for i in range(8):
        tr = int(r[i])
        if r[i]>tr+.5:
            rInt[i]=tr+1
            rDist[i] = (tr+1)-r[i]
        else:
            rInt[i]=tr
            rDist[i] = r[i]-tr
    if sum(rInt)%2==0:
        return rInt
    else:
        rDist.index(max(rDist))
        if rInt[i]>r[i]:
            rInt[i]=rInt[i]-1
        else: 
            rInt[i]=rInt[i]+1
        
    return rInt
    
def decodeE8(r):
    '''
        this decoder uses the D8 partition to decode E8 lattice
        it returns an integer label for the lattice point 
    '''
    nr = decodeD8(r)
    hr = [k+.5 for k in decodeD8([j-.5 for j in r])]
    
    nrd = sum([(nr[i]-r[i])**2 for i in range(8)])
    hrd = sum([(hr[i]-r[i])**2 for i in range(8)])
    
    #our range is 0-2, so normalized to 2
    if hrd>nrd:
        s = 0
        for i in range(8):
            s = s+((1<<(i*2)) *((nr[i]*2)%4))
        return s
    else:
        s = 0
        for i in range(8):
            s = s+((1<<(i*2)) *((hr[i]*2)%4))
        return s



randList = []

def decode24E8(v,randList):

    decoding = []#0
    for i in range(len(randList)/8):
        projVec = [0.0]*8
        for j in range(8):
            projVec[j]=(v[randList[i*8+j]]/8.0)*2.0
        decoding.append(decodeE8(projVec))# = decoding + (1<<(i*16))*decodeE8(projVec)

    return decoding
    



def testLatticeDecoders2():
    '''
        this test generates random vectors from a center point with increasing variance
    '''
    l=24
    
    rho=1.30
    repNumRnd=int(l**rho)+(8-int(l**rho)%8)
    randList1 = [randint(0,l-1) for i in range(repNumRnd)]
    
    rho=1.40
    repNumRnd=int(l**rho)+(8-int(l**rho)%8)
    randList2 = [randint(0,l-1) for i in range(repNumRnd)]
    
    rho=1.50
    repNumRnd=int(l**rho)+(8-int(l**rho)%8)
    randList3 = [randint(0,l-1) for i in range(repNumRnd)]

    rho=1.60
    repNumRnd=int(l**rho)+(8-int(l**rho)%8)
    randList4 = [randint(0,l-1) for i in range(repNumRnd)]

    import pylab
    numcntrs = 100
    numpts = numcntrs*1000
    maxvar = 1.5
    varlist = pylab.arange(0.01,maxvar,.05)
    #dec0s = [0.0]*len(varlist)
    dec1s = [0.0]*len(varlist)
    dec2s = [0.0]*len(varlist)
    dec3s = [0.0]*len(varlist)
    dec4s = [0.0]*len(varlist)
    
    
    
    
    
    avgdistanceList = [0.0]*len(varlist)
    cntrs = [[random()*8 for i in range(24)] for j in range(numcntrs)]
    i = 0
    
    for v in varlist:
        avgdistance = 0.0
        for m in range(numpts):
           ptcntr = int(random()*numcntrs)-1
           pt = [cntrs[ptcntr][j]+gauss(0.0,v) for j in range(24)]
           avgdistance =avgdistance +sum([(pt[j]-cntrs[ptcntr][j])**2 for j in range(24)])/float(numpts)
           #dec0s[i]= float(dec(pt)==dec(cntrs[ptcntr]))/float(numpts)   +dec0s[i]
           
           #special compare1
           counts = 0

           pDec = decode24E8(pt,randList1)
           cntrDec = decode24E8(cntrs[ptcntr],randList1)

           for j in range(len(pDec)):
                counts = int(pDec[j]==cntrDec[j])+counts
           
           dec1s[i]= int(counts>(2*len(pDec)/3.0))/float(numpts) +dec1s[i]
           
           
           
           
           #special compare2
           counts = 0
           pDec = decode24E8(pt,randList2)
           cntrDec = decode24E8(cntrs[ptcntr],randList2)
           for j in range(len(pDec)):
                counts = int(pDec[j]==cntrDec[j])+counts
           
           dec2s[i]= int(counts>(2*len(pDec)/3.0))/float(numpts) +dec2s[i]
           
           
           
           #special compare3
           counts = 0
           pDec = decode24E8(pt,randList3)
           cntrDec = decode24E8(cntrs[ptcntr],randList3)
           for j in range(len(pDec)):
                counts = int(pDec[j]==cntrDec[j])+counts
           
           dec3s[i]= int(counts>(2*len(pDec)/3.0))/float(numpts) +dec3s[i]
           
           
           
           
           
           #special compare4
           counts = 0
           pDec = decode24E8(pt,randList4)
           cntrDec = decode24E8(cntrs[ptcntr],randList4)
           for j in range(len(pDec)):
                counts = int(pDec[j]==cntrDec[j])+counts
           
           dec4s[i]= int(counts>(2*len(pDec)/3.0))/float(numpts) +dec4s[i]
           
        avgdistanceList[i] = avgdistance
        print v,avgdistanceList[i],dec1s[i],dec2s[i],dec3s[i],dec4s[i]
        i=i+1
       
    
    import pylab
    import math
    #pylab.plot(avgdistanceList,dec0s,label="Leech Lattice -.2dB")
    pylab.plot(avgdistanceList,dec1s,label="E8-Rand Proj-24 k=1.30")
    pylab.plot(avgdistanceList,dec2s,label="E8-Rand Proj-24 k=1.40")
    pylab.plot(avgdistanceList,dec3s,label="E8-Rand Proj-24 k=1.50")
    pylab.plot(avgdistanceList,dec4s,label="E8-Rand Proj-24 k=1.60")
    pylab.legend(loc='upper right')
    pylab.xlabel('Avgerage Distance')
    pylab.ylabel('Average Collisions Probability')

    pylab.show()




def genRandomCenters(numCntrs,rng=[0.0,8.0],length=24,spareseness=.90):
    '''
        this function will generate @numCntrs , random centers
        over the specified @rng range (defaul 0,7)
    '''
    centers = []
    for i in xrange(numCntrs):
        cntr = [0.0]*length
        
        for j in xrange(length):
            if random()>spareseness:
                cntr[j] = random()*rng[1]
        centers.append(cntr)
    return centers


def generateDataPts(numpts,centers=genRandomCenters(4),sigma = 0.20,length=24):
    '''
        this function generates data that is gaussian randomly perterbed from a random
        cluster center specified in centers.
        @numpts = takes a number of points to generate
        @centers = set of cluster centers
        @sigma = 4, variance of the normal variate
    '''
    data = [[] for i in range(numpts)]
    realCenters = [0]*numpts
    for i in xrange(numpts):
        pt = [0]*length
        #choose a cluster center that this point considers its representative
        clusterRep = randint(0,len(centers)-1)
        realCenters[i]=clusterRep
        for j in xrange(length):
            if centers[clusterRep][j]!=0.0:
                pt[j] = centers[clusterRep][j]+gauss(0,sigma)
        data[i]=pt
    return data,realCenters
    



testLatticeDecoders2()
#sim(10)













