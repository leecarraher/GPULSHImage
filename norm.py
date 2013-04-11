from pylab import rand,plot,show


def nn(p,q,d=2):
    ret = 0
    #for di in range(d):
    m = 10000.0
    for i in range(len(q)):
        d = (p[0]-q[i][0])**2+(p[1]-q[i][1])**2
        if d<m:
            m=d
            r=i   
    return r

def quant(n,r):
    l = []
    for i in range(1,n):
        x=r*2.0*i/n-r
        for j in range(1,n):
            y=2.0*r*j/n-r
            l.append((x,y))
    return l


def sim(numpts,numLat):
    print "started quantizer "+str(numLat)+"x"+str(numLat) 
    q = quant(numLat,1.0)
    
    #testpts
    R = (rand(numpts,2)*2.0)-1.0
    print "generated: "+str(numpts)+" random points"
    #(dist,count) - vector
    c = []
    for i in range(len(q)):
        c.append ((q[i][0]**2+q[i][1]**2,0 ))
    print "running simulation"
    #run sim
    
    for i in range(len(R)):
        l = nn(R[i],q)
        c[l]=(c[l][0],c[l][1]+1)
    c.sort()
    return c

        

c=sim(100000,4)
plot([cc[0] for cc in c],[cc[1] for cc in c])
c=sim(100000,8)
plot([cc[0] for cc in c],[cc[1] for cc in c])
c=sim(100000,12)
plot([cc[0] for cc in c],[cc[1] for cc in c])
c=sim(100000,16)
plot([cc[0] for cc in c],[cc[1] for cc in c])
show()
