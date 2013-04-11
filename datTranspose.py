import sys
infile = open(sys.argv[1],'r')
 
s = infile.readline()[:-1].split(" ")
l = len(s)
tran = [[]for h in range(l)]

while len(s)==l:
    for i in range(l):
        tran[i].append(float(s[i]))
    s = infile.readline()[:-1].split(" ")

if len(sys.argv)==3:
    outfile = open(sys.argv[2],'w')
    for i in range(l):
        outfile.write( str(tran[i]) + "\n")
else:
    for i in range(l):
        print str(tran[i])
