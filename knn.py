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


author: lee carraher

this is an ad-hoc version of the e2lsh (by Alexandr Andoni and 
Piotr Indyk) program in python it uses the leech decoder compiled
obj for the lsh function it then concatenates the hashes of the 
vectors for a given object at 240 bits and stores this with a 
pointer bach to the originating object in a hash table.

calling this module as main will build the hash table for items 
in the pgm directory and store it as storedDatabase.pickle

importing this module allows a stored database to be used for 
imaging matching unseen or test data, without having to 
regenerate the database each time, ie for 1 off tests or 
something like that. 

'''

import glob
keys75 = glob.glob("pgm/*.key")
import shelve


# leech library stuff
from ctypes import *
import os

import time,random

leechdec = cdll.LoadLibrary(os.getcwd()+'/leechDecoder.so')
#siftmatch = cdll.LoadLibrary(os.getcwd()+'/match.so')
#siftmatch.FindMatches2.argtypes= [POINTER(c_char),POINTER(c_char)]
#siftmatch.FindMatches2.restype = c_int

distC = cdll.LoadLibrary(os.getcwd()+'/dist.so')
distC.numMatches.argtypes= [c_int,c_int,POINTER(c_int),POINTER(c_int)]
distC.numMatches.restype = c_int

leechdec.decode.argtypes= [POINTER(c_float),POINTER(c_float)]
leechdec.decode.restype = c_int



floatArray =c_float*24

def decode(pts):
    cx=floatArray(*pts)
    d = c_float()
    result = leechdec.decode(cx,byref(d))

    return int(result),d.value
    
#def match(file1,file2):
#    f1=(c_char*len(file1))(*file1) 
#    f2=(c_char*len(file2))(*file2) 
#    result = siftmatch.FindMatches2(f1,f2)
#    return int(result)
    
    
def flatten(l):

    out = []
    for s in l:
        out.extend(s)
    return out

def distImgC(kfile1,kfile2):
    kfile1=kfile1
    kfile2=kfile2
    intArray1 =c_int*len(kfile1)
    intArray2 =c_int*len(kfile2)
    
    
    cx=intArray1(*kfile1)
    cy=intArray2(*kfile2)
    
    result = distC.numMatches(c_int(len(kfile1)),c_int(len(kfile2)),cx,cy)
    return int(result)
#end leech library stuff

l = 240
maxint = 0
def scale(d):
    for i in range(len(d)):
        d[i] = ((128-d[i]) /128.0)*8.0         

#shelves and dbs only hold strings and ints, we could naively 
#create strings of the longs, but then it assumes each is 8 bits
#this conversion is more efficient for storage space
def convert2Bytes(b):
    #return str(b)
    #print b[0]
    b = bin(b)[2:]
    ret = ''
    for i in range(len(b)/8):
        ret = ret + chr(int(b[i*8:(i+1)*8],2))
    return ret
    


import time

def decode128(keyvec,idxs,var=False):
    scale(keyvec)
    hashes = 0
    distance = 0.0
    starttime = 0
    tottime = 0
    
    for i in range(l/24):
        v = [0]*24
        for j in range(24*i,24*(i+1)):
            v[j%24] = keyvec[idxs[j]]
        starttime=time.time()
        dec = decode(v)
        tottime=tottime + (time.time()-starttime)
        distance = distance + 1.0/(dec[1]+1.0)
        hashes=hashes+ (dec[0] << 24*i)
    
    return hashes,distance,tottime

def dist(r,d):
    tot = 0.0
    for i in range(128):
        tot = tot+(r[i]-d[i])**2
    return tot**.5
    
    

def imageDist(img1list,img2list):
    '''
    returns the total distance between images
    '''
    totalDist =0.0
    if len(img2list)> 3000:### this is set because some images have >> than normal # of vectors
        return 100000.0
    
    for r in img1list:
        nearest = 1000000.0
        idx = 0
        for d in xrange(len(img2list)):

            l = dist(r,img2list[d])
            if l < nearest:
                nearest = l
                idx = d
            
        totalDist = totalDist + nearest
            
    return totalDist
    
#this should be iterative to avoid storing all of the keys in memory
def loadimages(files):
    data = []
    for fname in files:
        data.append((fname,loadimage(fname)))
    return data
    
def loadimage(fname):
    fileData = []
    f = open(fname,"r")
    f.readline()
    while f.readline() != "":
        vecData = []
        #should be 6 rows of 20 values, then one of 8 = 128
        for i in range(7):
            vecData.extend([int(k) for k in f.readline().split(" ")[1:]  ])
        fileData.append(vecData)
    f.close()
    return fileData

def hashcompare(v1,v2):
    incommon = 0
    for i in v1:

        for j in v2:
            if i ==j:incommon=incommon+1

    return incommon
    
    

def findExact(data):
    minTot = [0.0,""]
    for f in keys75:
        iD =  imageDist(data,loadimage(f))
        if iD < minTot[0] :
            minTot = [iD,f]
    print minTot
L=3
def findBestMatches(htable,data,idxs,k=10,exact=False):
    t = 0
    
    matches = {}
    for i in range(len(data[0])):
        
        for l in range(L):
            start = time.time()
            dec,distance,decTime = decode128(data[0][i],idxs,bool(l))
            
            #t = t+decTime
            vec = convert2Bytes(dec)
            t=t+(time.time()-start)
            if htable.has_key(vec):
                for pic in htable[vec]:
                    if matches.has_key(pic):
                        matches[pic] = matches[pic] +(distance)
                    else:
                        matches[pic] = (distance)

    i =0
    ret = [()]*len(matches.keys())        
    for key in matches.keys():
        ret[i] = (matches[key],key)
        i = i+1
        
    ret.sort(reverse=True)          
    
    '''
    matches = []
    # numvecs = len(data[0])
    for i in range(len(data[0])):
        dec = decode128(data[0][i],idxs)
        vec = convert2Bytes(dec[0])
        if htable.has_key(vec):matches.extend(htable[vec])
        
        
    s = set()
    for m in matches:
        s.add(m)
        
    ret = []
    #cnt = 0
    for p in s:
        #if cnt > k:break
        #cnt = cnt+1
        l = 1#len(loadimage(p))
        ret.append((matches.count(p),p))
        
    #for p in range(len(ret)):
    #    l = len(loadimage(ret[p][1]))
    #    ret[p] = (ret[p][0]/float(l),ret[p][1])
    ''' 
    
        
    #print time.time() -t
    #this is the additional step for exact comparison of candidates    
    if not exact:
        return ret[:k],t
        
    queryVectors = loadimage(data[1])
    exactDistances = []
    
    
    
    for distKeyURL in ret[0:k]:
        keyURL = distKeyURL[1]
        #print float(len(queryVectors))
        #print distKeyURL
        exactDistances.append((imageDist(queryVectors,loadimage(keyURL)), keyURL))
    
    exactDistances.sort()
    
    return exactDistances,t
        


import time
#in this case the htable is actually a shelf, persistent storage
def buildHashDatabase(files,htable,idxs):

    
    for fname in files:
        #read a sift vector file
        fileData = []
        f = open(fname,"r")
        f.readline()#just the file size
        while f.readline() != "":
            vecData = []
            #should be 6 rows of 20 values, then one of 8 = 128
            for i in range(7):
                vecData.extend([int(k) for k in f.readline().split(" ")[1:]  ])
            fileData.append(vecData)
        f.close()

        #hashing part
        fingerprint2 = []
        totalDecodeTime =0
        totalStartTime = 0
        for i in range(len(fileData)):#this is the conversion part
            startTime = time.time()
            decodeTime = time.time()
            vec = convert2Bytes(decode128(fileData[i],idxs)[0] )
            totalDecodeTime = totalDecodeTime + ( time.time() - decodeTime )
            if htable.has_key(vec):
                htable[vec].append(fname)
            else:
                htable[vec]=[fname]
            totalStartTime =totalStartTime+( time.time() -startTime)
        #print fname

        print  totalStartTime,totalDecodeTime
        
            
    print "done building database"
    return htable  
    
    
'''        
def buildHashDatabase(files,htable,idxs):


    for fname in files:
        #read a sift vector file
        fileData = []
        f = open(fname,"r")
        f.readline()#just the file size
        while f.readline() != "":
            vecData = []
            #should be 6 rows of 20 values, then one of 8 = 128
            for i in range(7):
                vecData.extend([int(k) for k in f.readline().split(" ")[1:]  ])
            fileData.append(vecData)
        f.close()

        #hashing part
        fingerprint2 = []
        for i in range(len(fileData)):
            vec = decode128(fileData[i],idxs) 
            if htable.has_key(vec):
                htable[vec].append(fname)
            else:
                htable[vec]=[fname]
        print fname
            
    print "done building database"
    return htable
'''

if __name__ == '__main__':
    pass

    import random
    #analysis shows that the prob of selecting all vectors at least once is 
    #(1-(1-1/m)^n)^m which for n=240 and m=128 is very small 6.6e-10
    #480 -> .05 , 672-> .5 and 960 is 95%, but these vectors are too long
    #so we will force it here
    idxs = range(128)
    idxs.extend([random.randint(0,127) for i in range(128,l)])
    random.shuffle(idxs)


    #htable = {}#
    htable = shelve.open("datTable.dat","n",writeback=True)
    totaltime = time.time()
    buildHashDatabase(keys75,htable,idxs)

    print "total time: "+str(time.time()-totaltime)
    
    #import gzip
    import cPickle
    cPickle.dump(idxs,file("idxs.dat","w"))
    #f.close()

