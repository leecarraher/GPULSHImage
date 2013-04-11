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

import time,random,numpy,math,os

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
l = 2
tblSize = 4096**l



def cnt(v):
    c =  ((v & 0xfff) * 0x1001001001001 & 0x84210842108421) % 0x1f
    c += (((v & 0xfff000) >> 12) * 0x1001001001001 & 0x84210842108421)%0x1f
    return c
    
maxlen =0
'''
maybe a cleverway to partition hashes
using what we know about octads
2576^n options
2*759^n
2*2^n
and all permutations
'''    
'''def ELFHash(key):
    #key=str(key)
    hash = 0
    x = 0
    while key>0:
        part = key&((1<<24)-1)
        key = key>>24
        
        t = cnt(part)
        #largest subset is 12
        if t==12:
            hash = (hash<<12)+(part%4079)
        elif t ==0 or t ==24:
            hash = (hash<<10)+(part%751)
        else:
            hash = (hash<<1)+(part%2)

    return hash
'''

def ELFHash(key):
    dist = key[1] 
    key = key[0]
    key = str(key)
    hash = 0
    x    = 0
    for i in range(len(key)):
      hash = (hash << 4) + ord(key[i])
      x = hash & 0xF0000000
      if x != 0:
        hash ^= (x >> 24)
      hash &= ~x
    return hash%tblSize,dist





def scale(d,low,high,newLow,newHigh):
    k = (newHigh - newLow) / (high - low)
    for i in xrange(len(d)):        
        d[i] = newLow+(d[i]-low)*k



def decode(pts):

    scale(pts,-1.0,1.0,-8.0,8.0)

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
#constant scaling factor 1/sqrt(24)
sc = 1.0/((24)**.5)
maxint = 0

#shelves and dbs only hold strings and ints, we could naively 
#create strings of the longs, but then it assumes each is 8 bits
#this conversion is more efficient for storage space
def convert2Bytes(b):
    b = bin(b)[2:]
    ret = ''
    for i in range(len(b)/8):
        ret = ret + chr(int(b[i*8:(i+1)*8],2))
    return int(ret)


def decodegt24(keyvec,m,var=False):
    hashes = 0
    distance = 0.0
    #print sum(keyvec)/len(keyvec)
    scale(keyvec,0.0,170.0,-1.0,1.0)
    #print sum(keyvec)/len(keyvec)
    for i in xrange(len(m)):
        
        if(var):v = [k+random.gauss(0,1)/sc for k in numpy.dot(keyvec,m[i]).tolist() ]
        else:v = numpy.dot(keyvec,m[i]).tolist()
        dec = decode(v)
        distance = distance + 1.0/(dec[1]+1.0)
        hashes=hashes+ (dec[0] << 24*i)
    return hashes,distance

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

def query(htables,data,mats,k=10,exact=False,sort=True):
    t = 0
    rho = k**(.3641)#len(htables[0])**(.3641)#from monte carlo simulations
    matches = {}
    
    #iterate over the vectors in an image
    for i in xrange(len(data)):
        for j in xrange(len(htables)):
            for l in xrange(int(rho)):
                dec,distance = ELFHash(decodegt24(data[i],mats[j]) )
                if not htables[j][dec]==None:
                    for pic in htables[j][dec][0:10]:
                        if matches.has_key(pic):
                            matches[pic] = min(matches[pic] ,(distance)) #(distance+matches[pic])/2.0#
                        else:
                            matches[pic] = (distance)

    i =0
    ret = [()]*len(matches.keys())        
    for key in matches.keys():
        ret[i] = (matches[key],key)
        i = i+1
        
    if sort:ret.sort()          

        
    #print time.time() -t
    #this is the additional step for exact comparison of candidates

    if not exact:
        return ret
        
    queryVectors = loadimage(data[1])
    exactDistances = []
    
    
    
    for distKeyURL in ret[0:k]:
        keyURL = distKeyURL[1]
        #print float(len(queryVectors))
        #print distKeyURL
        exactDistances.append((imageDist(queryVectors,loadimage(keyURL)), keyURL))
    
    exactDistances.sort()
    
    return exactDistances,t
        

#in this case the htable is actually a shelf, persistent storage
def buildHashDatabase(files,htable,mats):
    maxlen =0
    inc = 0
    fcount =0
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
        for i in range(len(fileData)):
            vec = ELFHash(decodegt24(fileData[i],mats) )[0]
            if htable[vec]==None:
                inc=inc+1
                htable[vec]=[fcount]
            else:
                htable[vec].append(fcount)

                

        fcount = fcount+1
        #print fname,inc,fcount
    
    return htable  
    
def buildHashDatabase2(files,htable,mats):
    maxlen =0
    inc = 0
    fcount =0
    for fileData in files:



        #hashing part
        fingerprint2 = []
        totalDecodeTime =0
        totalStartTime = 0
        for i in range(len(fileData)):
            vec = ELFHash(decodegt24(fileData[i],mats) )[0]
            if htable[vec]==None:
                inc=inc+1
                htable[vec]=[fcount]
            else:
                htable[vec].append(fcount)

                

        fcount = fcount+1
        #print fname,inc,fcount
    
    return htable  
        

if __name__ == '__main__':
    pass
    import cPickle
    repeats=5
    

    for m in range(repeats):
        htable =[ None  for i in xrange(tblSize)]#shelve.open("datTableMat"+str(m)+".dat","n")
        mats = [numpy.array([[random.gauss(0,1)*sc for i in range(24)]for j in range(128)]) for k in range(l)]
        totaltime = time.time()
        buildHashDatabase(keys75,htable,mats)
        cPickle.dump(mats,file("randomMatD2"+str(m)+".dat","w"))
        cPickle.dump(htable,file("datTableMatD2"+str(m)+".dat","w"))
        print "done building database " +str(m+1)
    print "total time: "+str(time.time()-totaltime)
    


