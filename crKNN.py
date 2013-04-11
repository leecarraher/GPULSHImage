'''

Copyright 2010 Lee Carraher. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Lee Carraher ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Lee Carraher OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of Lee Carraher.


a very simple test module that runs all the dirtied image siftkeys agains the clean sift keys
outputs results of searches, by finding the k=20 nearest neighbors
#gcc -c -fPIC dist.c
#gcc -shared dist.o -o dist.so

'''
# leech library stuff
from ctypes import *
from math import log,sqrt
from numpy import array,dot
import os
from random import gauss,random,randint

leechdec = cdll.LoadLibrary(os.getcwd()+'/leechDecoder.so')
leechdec.decode.argtypes= [POINTER(c_float),POINTER(c_float)]
leechdec.decode.restype = c_ulonglong

floatArray =c_float*24


def scale(d,low,high,newLow,newHigh):
    k = (newHigh - newLow) / (high - low)
    for i in xrange(len(d)):        
        d[i] = newLow+(d[i]-low)*k

def distance(x,y,d):
    dist =0.0
    for i in xrange(d):
        dist = dist+(x[i]-y[i])*(x[i]-y[i])
    return dist**.5
    

def linearSearchNN(q,DB,d,minVal=0.0,maxVal=1.0):
    minDist = 1000000.0
    
    for i in xrange(len(DB)):
        dist = distance(q,DB[i],d)
        if dist<minDist:
            minDist = dist
            minq = i
    return minq



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
    scale(r,0.,1.,0.0,2.0)

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
        
def decode(pts):
    scale(pts,-1.0,1.0,-8.0, 7.0)
    cx=floatArray(*pts)
    d = c_float()
    result = leechdec.decode(cx,byref(d))
    #print bin(result)
    return result,d.value
    
#end leech library stuff

#constant scaling factor 1/sqrt(24)
sc = 1.0/(sqrt(24))
maxint = 0


def convert2Bytes(b):
    '''
    shelves and dbs only hold strings and ints, we could naively 
    create strings of the longs, but then it assumes each is 8 bits
    this conversion is more efficient for storage space
    '''
    b = bin(b)[2:]
    ret = ''
    for i in range(len(b)/8):
        ret = ret + chr(int(b[i*8:(i+1)*8],2))
    return ret
    
    
    
def generateRandomProj(k,sim,d):
    M = [[] for i in xrange(sim)]
    for i in range(sim):
        for j in range(k):
            M[i].append( array([[gauss(0,1)*sc for a in range(24)] for b in range(d)]))
    return M
    

def decodeGt24(v,M,minVal,maxVal,rand=False):
    scale(v,minVal,maxVal,0.0,1.0)
    distance = 0.0
    hashVal = 0
    
    for i in xrange(len(M)):
        vproj = dot(v,M[i]).tolist()
        if rand:vproj=[vhat+gauss(0,.047)/sc for vhat in vproj]
        dec = decode(vproj)

        distance = distance + 1.0/(dec[1]+1.0)
        hashVal=dec[0]#hashVal+ (dec[0]<< (24*i))
    return convert2Bytes(hashVal)


def buildDB(data,minVal=0.0,maxVal=1.0):
    n = len(data)
    d = len(data[0])
    #parameters
    P1=.01779
    P2=.0000156
    rho = log(P1)/log(P2)
    k =  int((log(n)/log( 1/P2)+.5))
    sim = int(log(n)+.5)
    
    M = generateRandomProj(k,sim,d)
    
    print n**rho,k,sim
    
    DBs = [{} for s in xrange(sim)]
    for s in xrange(sim):
        for i in xrange(n):
            hashVal = decodeGt24(data[i],M[s],minVal,maxVal)
            if DBs[s].has_key(hashVal):
                DBs[s][hashVal].append(i)
            else:
                DBs[s][hashVal]=[i]
                    
    return DBs,M
    
def query(q,DBs,Ms,n,l=1000,minVal=0.0,maxVal=1.0):
    #parameters
    P1=.01779
    P2=.0000156
    rho = log(P1)/log(P2)
    sims = len(Ms)#number of multi-runs for whp
    k = len(Ms[0])#number of random matrix projections per run
    
    candidates = set()
    #first iterate over the runs
    for s in xrange(sims):
        #next iterate over the n^rho nearby points
        hashVal = decodeGt24(q,Ms[s],minVal,maxVal)
        if DBs[s].has_key(hashVal):
        
            for c in DBs[s][hashVal]:
                candidates.add(c)
                
        for r in xrange(int(n**rho+.5)):
            hashVal = decodeGt24(q,Ms[s],minVal,maxVal,True)
            
            if DBs[s].has_key(hashVal):
                for c in DBs[s][hashVal]:
                    candidates.add(c)
            if len(candidates)>2*l:return candidates
            
    return candidates





#<-----------Testing stuff--------------->    
import pylab
def testCollisionPr(n ,d=24):
     
    M =array([[gauss(0,1)/(24**.5) for a in range(24)] for b in range(d)])
    S = [0.0]*n
    C = [0]*n
    #generate distances and buckets
    i=0
    while i <n:

        p = [(random()*2.0)-1.0 for j in xrange(d)]
        q = [p[j] + (gauss(0.0,1.0)/(d**.5)) for j in xrange(d)]
        S[i]=distance(p,q,d)
        C[i]= decode(p)[0]==decode(q)[0]  
        i = i+1
        
    ranges = pylab.histogram(S,40)[1]   
    bucketsCol = [0]*len(ranges)
    bucketsDis = [0]*len(ranges)


    #fill buckets with counts 
    for i in xrange(n):
        k = len(ranges)-1
        while S[i] < ranges[k]:k=k-1
        if C[i]:bucketsCol[k]=bucketsCol[k]+1
        bucketsDis[k] = bucketsDis[k]+1
        
    print bucketsDis
    print bucketsCol
    print ranges

    p = [0]*len(ranges)
    for m in range(len(ranges)):
        if bucketsDis[m]>0:
            p[m]=bucketsCol[m]/float(bucketsDis[m])

    pylab.plot(ranges,p) 


def testCollisionsE8(n,d=24):
    M = eye(24,24)

    S = [0.0]*n
    C = [0]*n
    #generate distances and buckets
    for i in range(n):
        p = [random() for j in xrange(d)]
        q = [p[j] + (gauss(0,1)/(d**.5)) for j in xrange(d)]
        S[i]=distance(p,q,d)
        C[i]= int(decodeE8(dot(p,M)) == decodeE8(dot(q,M)))
    
    ranges = pylab.histogram(S,30)[1]   
    bucketsCol = [0]*len(ranges)
    bucketsDis = [0]*len(ranges)

    #fill buckets with counts 
    for i in xrange(n):
        k = len(ranges)-1
        while S[i] < ranges[k]:k=k-1
        if C[i]:bucketsCol[k]=bucketsCol[k]+1
        else:bucketsDis[k] = bucketsDis[k]+1
    print bucketsDis
    print ranges
    pylab.plot(ranges,[float(bucketsCol[i])/(float(bucketsDis[i]+.000000000001))  for i in range(len(ranges))]) 
    
    
def compareToLinearAcc(n=1000):
    d=120
    dsc = 1/sqrt(d)
    
    randVectors = [[random() for j in xrange(d)] for i in xrange(n)]
    DBs,Ms = buildDB(randVectors)
    print [len(DBs[i]) for i in  range(len(DBs))]
    
    counts = []
    testLen = 1000
    distances = 0.0
    for i in xrange(testLen):
        q = randint(0,n-1)
        candidate = randVectors[q]
        realNN = linearSearchNN(randQuery,randVectors,d)
        approxNNs = query(randQuery,DBs,Ms,n)
        for approxNN in approxNNs:
            if realNN==approxNN:   
                print len(approxNNs)
       
    import pylab
    gram = pylab.histogram(counts)
    pylab.plot(gram[1][:-1],gram[0]) 




if __name__ == '__main__':
    pass
    #testCollisionsE8(400000)
    testCollisionPr(8000000)
    pylab.show()

    
