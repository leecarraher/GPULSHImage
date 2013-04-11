import cPickle as pickle
import shelve
import time
start = time.time()
idx,htable = pickle.load(open("storedDatabase.pickle","r"))
print "done loading pickle: "+ str(time.time()-start)

#shelves and dbs only hold strings and ints, we could naively 
#create strings of the longs, but then it assumes each is 8 bits
#this conversion is more efficient for storage space
def convert2Bytes(b):
    #return str(b)
    b = bin(b)[2:]
    ret = ''
    for i in range(len(b)/8):
        ret = ret + chr(int(b[i*8:(i+1)*8],2))
    return ret

#unique filenames
    '''s = set()
    for key in htable.keys():
        for f in htable[key]:
            s.add(f)
     
    print "filename set creates"       
    #filename->unique#
    name2num = {}
    i = 0
    for name in s:
        name2num[name]=i
        i=i+1
    print "filenames mapped to unique ints" 
''' 

newHtable = shelve.open("datTable.dat","n",writeback=True)
#newHtable = {}
for key in htable.keys():
    newHtable[convert2Bytes(key)]=htable[key]
    '''
        fulllist = htable[key]
        t = []
        for name in fulllist:
            t.append(name2num[name])
        newHtable[convert2Bytes(key)] = t
    print "newHtable created, filenames replaced by unique nums"
        
    del(htable)

    namelst = ['']*len(s)
    for name in s:
        namelst[name2num[name]] = name
    '''
print "dumping the pickle"
pickle.dump(idx,open("idxs.dat",'w'))
newHtable.sync()
newHtable.close()

#start =time.time()
#idx,namelst,htable = pickle.load(open("anotherpickle.pck","r"))
#print "done loading pickle: "+str(time.time()-start)
    
