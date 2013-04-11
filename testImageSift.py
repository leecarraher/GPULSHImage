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
from knn import *
import cPickle as pickle
print "loading image files"
import glob
keys15 = glob.glob("pgmdirty/*.key")


    
#really slow, needs to be fixed
#dat15=loadimages(keys15)

#import gzip
import shelve


print "loading database"
#f = open('storedDatabase.pickle','r')
#idxs,htable = pickle.load(f)
idxs = pickle.load(file("idxs.dat","r"))
htable = shelve.open("datTable.dat","r")

#f = gzip.open('storedDatabase.pickle.gz','rb')
#idxs,htable = pickle.load(f)

print len(htable)
totaltime = 0.0
starttime = -1.0
totaltimeDec = 0.0
import time

correct20 = 0
correct10 = 0
correct5 = 0
correct = 0
correct0 = 0

def loadAllKeys(keys1,keys2):
    
    allKeys75 = []
    allKeys15 = []
    for i in xrange(len(keys1)):
        l = loadimage(keys1[i])
        m = loadimage(keys2[i])
        if len(l)<400 and len(m):
            allKeys15.append(flatten(l))
            allKeys75.append(flatten(m))
    return allKeys15,allKeys75



klist15,klist75 = loadAllKeys(keys15,keys75)

#import pylab
#graph = []
#for k in klist:
#    graph.append(len(k))

#graph.sort()
#pylab.plot(graph)
#pylab.show()

print len(klist15)
print len(klist75)

for i in range(0,len(keys15)):
    starttime = time.time()
    data = (loadimage(keys15[i]),keys15[i])
    
    starttime = time.time()
    #argmin
    mostMatches = [0,0]
    mname = keys15[i].split("/")[1]
    
    query = klist15[i] 
    for j in xrange(len(klist75)):
        m = distImgC(query,klist75[j])
        #if j%50==0:print j
        if m>mostMatches[0]:
            mostMatches[0]=m
            mostMatches[1]=j
            
    tempName = keys75[mostMatches[1]].split("/")[1]
    mname = keys15[i].split("/")[1]
    if mname==tempName:
        correct0 = correct0 +1
    #print "\t realmatches: "+ str(distImgC(query,klist75[i]))  + " file: "+keys75[i]
    #print "\t argmin All: " + str(mostMatches[0]) + " file: "+keys75[mostMatches[1]]
    
    ltime= time.time() - starttime
    
    
    
    starttime = time.time()
    # there is an option for exact matching, but it takes an extremely long time to check just a few files
    t,timeDec = findBestMatches(htable,data,idxs,20)
    #print timeDec,time.time()-starttime
    #totaltime = totaltime + (time.time()-starttime)
    #totaltimeDec = totaltimeDec + timeDec
    
    #objtype = mname.split("_")[0]
    #types = 0
    for j in range(len(t)):
        tempName = t[j][1].split("/")[1]

        if tempName==mname:correct20 = correct20+1
        if tempName==mname and j <=15:correct10 = correct10+1
        if tempName==mname and j<=10:correct5 = correct5+1
        if tempName==mname and j<=1:correct = correct+1
        #if tempName.split("_")[0]==objtype:types=types+1
        
    print i,correct20,correct10,correct5,correct,correct0,time.time()-starttime,ltime
    

print "total images = " + str(len(keys15))
print "total time = " + str(totaltime)
print "total decoding time = " + str(totaltimeDec)
print "Average Time = " + str(totaltime/float(len(keys15)))
print "top20\ttop10\ttop5\t1st"
print correct20,correct10,correct5,correct

