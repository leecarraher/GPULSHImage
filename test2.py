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
from ctypes import *
distC = cdll.LoadLibrary(os.getcwd()+'/dist.so')
distC.numMatches.argtypes= [c_int,c_int,POINTER(c_int),POINTER(c_int)]
distC.numMatches.restype = c_int

def distImgC(kfile1,kfile2):
    intArray1 =c_int*len(kfile1)
    intArray2 =c_int*len(kfile2)
    
    
    cx=intArray1(*kfile1)
    cy=intArray2(*kfile2)
    
    result = distC.numMatches(c_int(len(kfile1)),c_int(len(kfile2)),cx,cy)
    return int(result)

'''
from knn2 import *
import cPickle as pickle
print "loading image files"
import glob
keys15 = glob.glob("pgm/*.key")


    
#really slow, needs to be fixed
#dat15=loadimages(keys15)

#import gzip
#import shelve




totaltime = 0.0
starttime = -1.0
totaltimeDec = 0.0
import time

correct20 = 0
correct10 = 0
correct5 = 0
correct = 0
correct0 = 0
NScorrect0 = 0

def loadAllKeys(keys1,keys2):
    
    allKeys75 = []
    allKeys15 = []
    trunKeys = []
    for i in xrange(len(keys1)):
        l = loadimage(keys1[i])
        m = loadimage(keys2[i])
        if len(l)<400 and len(m):
            allKeys15.append(flatten(l))
            allKeys75.append(flatten(m))
            trunKeys.append(keys1[i])
    return allKeys15,allKeys75,trunKeys


def SelectTest():
    print "loading database"
    htables = []
    mats = []
    for m in glob.glob("randomMatD2*"):
        mats.append( pickle.load(file(m)))
    for h in glob.glob("datTableMatD2*"):
        htables.append( pickle.load(file(h)))
        #klist15,klist75,trunKeys = loadAllKeys(keys15,keys75)
    fulltest = False


    ctall = 0

    for name in keys15: 
        mname = name.split("/")[1]
        data = loadimage(name) 
        #approx Knn-LSH
        starttime = time.time()
        # there is an option for exact matching, but it takes an extremely long time to check just a few files
        t = query(htables,data,mats)
        
        tim = time.time()-starttime
        location = -1
        for j in xrange(len(t)):
            tempName = keys15[t[j][1]].split("/")[1]
            if tempName==mname and location==-1:
                location = j
                correct0 = correct0+1
        
        ctall = ctall +1
        
        
        
        print ctall, correct0,len(t) , location,tim#, len(NSt),NSloc,NStim
        
        #print len(t),correct20,correct10,correct5,correct,correct0,loc,tim
         

    
def online_variance(data,n=0,mean = 0,M2=0):

    for x in data:
        n = n + 1
        delta = x - mean
        mean = mean + delta/float(n)
        M2 = M2 + delta*(x - mean) 
        
    if n ==0:
        print "error"
        return 0,n,mean,M2
    
    variance_n = M2/float(n)
    variance = M2/float((n - 1))
    return variance,n,mean,M2
import random
import copy
def timeTest(start,stop,step): 
    print " N  \tLinear Mean  \t LSH Mean \tLinear Variance \tLSH Variance " 
    for size in range(start,stop,step):

        fileset = []
        lSearchDB = []
        for j in range(max(size,len(keys15))):
            name = keys15[j]
            l = loadimage(name)
            lSearchDB.append(l)
        size = len(lSearchDB)
        ltimes = []
        times = []
        #we have to fix the query DB here, its way bigger than the linear search db
        newHtables = []
        newMats = []
        l = 2

        for m in range(2):
            htable =[ None  for i in xrange(tblSize)]#shelve.open("datTableMat"+str(m)+".dat","n")
            mats = [numpy.array([[random.gauss(0,1)*sc for i in range(24)]for j in range(128)]) for k in range(l)]
            buildHashDatabase2(copy.deepcopy(lSearchDB),htable,mats)
            newHtables.append(htable)
            newMats.append(mats)
        #done fixing
        randQs = range(size)
        random.shuffle(randQs)
        randQs = randQs[:50]
        #select just a few (50) to test
        for n in randQs:
            name = keys15[n]
            mname = name.split("/")[1]
            data = lSearchDB[n]
            
            #exact search
            mostMatches = [0,0]
            starttime = time.time()
            '''
            for j in xrange(size):
                
                m = distImgC(flatten(data),flatten(lSearchDB[j]))
                if m>mostMatches[0]:
                    mostMatches[0]=m
                    mostMatches[1]=j
            '''              
            ltimes.append( time.time() - starttime)
            #approx Knn-LSH
            starttime = time.time()
            t = query(newHtables,copy.deepcopy(data),newMats,size)
            times.append(time.time()-starttime)
            
        
        variance,n,mean,M2 = online_variance(ltimes)
        
        lshvariance,n,lshmean,M2 = online_variance(times)
        print  str(size)+"\t"+str(mean)+ "\t"+str(lshmean)+ "\t"+str(variance)+"\t"+str(lshvariance)


    
#SelectTest()
timeTest(900,3000,100)


