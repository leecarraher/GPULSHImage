'''Copyright 2010 Lee Carraher. All rights reserved.

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
'''
'''
binary golay 24 12 8, encoder and imld and nearest neighbor decoder (fixed all working)
'''
I = [[1,0,0,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0,0,0 ],[0,0,0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,0,0,1]]
H = [[0,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,0,1,1,1,0,0,0,1,0],[1,1,0,1,1,1,0,0,0,1,0,1],[1,0,1,1,1,0,0,0,1,0,1,1],[1,1,1,1,0,0,0,1,0,1,1,0],[1,1,1,0,0,0,1,0,1,1,0,1],[1,1,0,0,0,1,0,1,1,0,1,1],[1,0,0,0,1,0,1,1,0,1,1,1],[1,0,0,1,0,1,1,0,1,1,1,0],[1,0,1,0,1,1,0,1,1,1,0,0],[1,1,0,1,1,0,1,1,1,0,0,0],[1,0,1,1,0,1,1,1,0,0,0,1]]
Hdec = [4094,1143,2619,3357,1679,2887,3491,3793,1897,949,475,2285]
Idec = [1,2,4,8,16,32,64,128,256,512,1024,2048]
words = []

#dec to binary
def D2B(n):
	bStr = [0]*12
	i = 0
	while n > 0:
		bStr[i]=n % 2
		n = n >> 1
		i = i+1
	return bStr

#binary to dec
def B2D(n):
	d = 0
	s =1
	for i in n:
		if i ==1:
			d=d+s
		s=s<<1
	return d

def encode(s):
	c = s
	#first part is just the identity
	c.extend([0]*12)
	#parity matrix
	for i in range(12):
		sum = 0.0
		for j in range(12):
			sum=sum+H[j][i]*s[j]
		c[i+12] = int(sum%2)	
	return c

#nearest neighbor
def nn(c):
	minimum = 24
	s = [0]*24
	for w in words:
		dist =0
		for i in range(24):
			if w[i]!=c[i]:
				dist=dist+1
		if dist<minimum:
			minimum=dist
			s=w
	return s

def add(c,u):
	l = [0]*len(c)
	if len(c) != len(u):
		print "ERROR!"
		return []
	for i in range(len(c)):
		l[i] = (c[i]+u[i])%2
	return l

def decode(c):
	'''
	Based on IMLD in Dekker Coding Theory: The Essentials
	'''
	#compute the 1st syndrome
	s1 = [0]*12
	for i in range(12):
		su = 0
		for j in range(12):
			su = su+I[i][j]*c[j]
		for j in range(12):
			su=su+H[i][j]*c[j+12]
		s1[i]=int(su%2)

	wts1 = sum(s1)
	if wts1 == 0:return c
	if wts1 <=3:
		s1.extend([0]*12)
		return add(c,s1)

	#2nd test for wt >3
	for i in range(12):
		s1H = add(s1,H[i])
		wt = sum(s1H)
		if wt <=2:
			u = [0]*24
			u[i] = 1
			for j in range(12):u[j]=s1H[j]%2
			return add(c,u)
	#2nd syndrome
	s2 = [0]*12
	for i in range(12):
		su = 0
		for j in range(12):su = su+H[i][j]*s1[j]
		s2[i]=int(su%2)
	
	wts2 = sum(s2)

	if wts2 <=3:
		u = [0]*12
		u.extend(s2)
		return add(c,u)
	
	#2nd test for second syndrom wt >3
	for i in range(12):
		s2H = add(s2,H[i])
		wt = sum(s2H)
		tu = [0]*24
		tu[i] = 1
		for j in range(12):
			tu[j+12]=s2H[j]%2

		if wt <=2:
			u = [0]*24
			u[i] = 1
			for j in range(12):
				u[j+12]=s2H[j]%2
			return add(c,u)

from random import randint
#build syndrome table
for i in range(0,2**12):
	words.append(encode(D2B(i)))

for i in range(20):
	w = randint(0,12)
	c = [0]*12
	for j in range(w):
		c[randint(0,11)] = 1

	print ""
	print "OW  " + str(c)
	s = encode(c)
	print "EW  " + str(s)

	r = randint(0,23)
	s[r] = (s[r]+1)%2 # flip this bit

	r = randint(0,23)
	s[r] = (s[r]+1)%2 # flip this bit

	r = randint(0,23)
	s[r] = (s[r]+1)%2 # flip this bit


	print "RW  " + str(s)
	print "DW  " + str(decode(s))
	print "nnDW" + str(nn(s))
