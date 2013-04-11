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
'''
from knn import *
import os
import cPickle as pickle
from tkFileDialog import askdirectory
from Tkinter import *
from PIL import *
import ImageTk
import Image
import shelve

print "loading database"
#f = open('storedDatabase.pickle','r')
#idxs,htable = pickle.load(f)
idxs = pickle.load(file("idxs.dat","r"))
htable = shelve.open("datTable.dat","r")


tempFiles = []

#store this for later comparison
compareKeyfile = ""
# uses sift and grabs the vectors 
def getSifts(src):
    if src.endswith("jpg"):
        os.system("convert"+ " " + src + " " + src + ".pgm")
        src = src+".pgm"
        tempFiles.append(src)

    print src,src+".key"
    os.system("./sift"+ " <" + src + "> " + src + ".key")
    tempFiles.append(src + ".key")
    vs = loadimage(src + ".key")
    global compareKeyfile
    compareKeyfile = src + ".key"

    return vs

#the os call to ./match for the queryImg and matchedImg
def createCompareImg(queryImg,matchedImg):
        if matchedImg.endswith("jpg"):
            matchedImg = matchedImg+".pgm"
            
        #print "./match -im1 "+ queryImg+" -k1 "+compareKeyfile+" -im2 "+matchedImg+" -k2 "+matchedImg+".key > out.pgm"
        os.system("./match -im1 "+ queryImg+" -k1 "+queryImg+".key"+" -im2 "+matchedImg+" -k2 "+matchedImg+".key > out.pgm")
        
    

root =Tk()
root.title('Find Similar Images')
filename = StringVar()
filters = StringVar()


#generate sift vectors for the query image and query the database for those vectors
def search(filename):

    vs = getSifts(filename)
    return findBestMatches(htable,(vs,filename),idxs,24)
    


previmg = 0
imgList = []

#load 5 images from the current directory
def loadimgs(opt=1):
    global imgList
    #previmg = 0
    imgList = glob.glob(filename.get()+'/'+filters.get()+'*.jpg')
    imgList.extend(glob.glob (filename.get()+'/'+filters.get()+'*.pgm') )
    
    for i in range(previmg,previmg+5):
        if len(imgList)>i:
            b = createButton(imgList[i],f2,i)
    
    
#query the match database and load candidate images
queryImage = ""
matches = []
def loadMatches(num):
    for img in f3.children.keys():
        f3.children[img].destroy()
    global queryImage
    global matches
    queryImage = imgList[num]
    
    matches = search(imgList[num])
    for num in range(len(matches)):
        filename = matches[num]
        
        filename = filename[1][:-4]#remove the .key part
        
        createButtonMatches(filename,f3,num,matches[num][0])

#program invokation of ./match for generating the complete vector match image for 2 images
def compareFiles(num):
    createCompareImg(matches[num][1][:-4],queryImage)
    im = Image.open("out.pgm")
    im.show()

#generate image buttons for 5 images in the query image directory
def createButton(filename,frame,num,w=150,h=75):
    
    image = Image.open(filename)
    image = image.resize((w, h),Image.ANTIALIAS)
    photo1 = ImageTk.PhotoImage(image=image)
    button1 = Button(frame, compound=TOP, width=w, height=h+7, image=photo1,
        text=filename.split('/')[-1][0:12], command=lambda: loadMatches(num))
    button1.grid(column=num%5+1 ,row=0)
    button1.image = photo1

#generate an image button for a candidate image, placed in a 6xn grid
def createButtonMatches(filename,frame,num,sims,w=150,h=75):
    image = Image.open(filename)
    image = image.resize((w, h),Image.ANTIALIAS)
    photo1 = ImageTk.PhotoImage(image=image)
    button1 = Button(frame, compound=TOP, width=w, height=h+7, image=photo1,text=str(sims)+":"+filename.split('/')[-1][0:12], command=lambda: compareFiles(num))
    button1.grid(column=num%6,row=int(num/6.0))
    button1.image = photo1
    
    
#function binding for displaying the next 5 images in the directory
def nextImgs():
    for img in f2.children.keys():
        f2.children[img].destroy()
    global previmg
    previmg = previmg+5
    loadimgs()
    btnnextImgs = Button(f2, text='Next >>',command=nextImgs)
    btnnextImgs.grid(column=6,row=0,padx=2,pady=2)
    btnprevImgs = Button(f2, text='<< Prev',command=prevImgs)
    btnprevImgs.grid(column=0,row=0,padx=2,pady=2)
    
#function binding for displaying the previous 5 images in the directory
def prevImgs():
    for img in f2.children.keys():
        f2.children[img].destroy()
        
    global previmg
    previmg = previmg-5
    loadimgs()
    btnnextImgs = Button(f2, text='Next >>',command=nextImgs)
    btnnextImgs.grid(column=6,row=0,padx=2,pady=2)
    btnprevImgs = Button(f2, text='<< Prev',command=prevImgs)
    btnprevImgs.grid(column=0,row=0,padx=2,pady=2)

#pop up a directory browser
def browser():
    filename.set(askdirectory())
    global previmg
    previmg = 0
    loadimgs()


#this is the directory selection frame
f1 = LabelFrame(root, text="Select Image Directory",height=90,width=150*7)
f1.pack(fill="both", expand="yes")

e = Entry(f1, width=60,textvariable=filename)
e.pack(side=LEFT,padx=2,pady=2)

e.bind('<Return>', loadimgs)


e2 = Entry(f1, width=10,textvariable=filters)
e2.pack(side=RIGHT,padx=2,pady=2)
e2.bind('<Return>', loadimgs)
e2Txt = Label(f1,text="Filter Images")
e2Txt.pack(side=RIGHT,padx=2,pady=2)

openDir = Button(f1, text='Load...',command=loadimgs)
openDir.pack(side= LEFT,padx=2,pady=2)
browse = Button(f1, text='Browse',command=browser)
browse.pack(side= LEFT,padx=2,pady=2)

#this is the query image fram
f2 = LabelFrame(root, text="Images",height=90,width=150*7)
f2.pack(fill="both", expand="yes")

btnnextImgs = Button(f2, text='Next >>',command=nextImgs)
btnnextImgs.grid(column=6,row=0,padx=2,pady=2)

btnprevImgs = Button(f2, text='<< Prev',command=prevImgs)
btnprevImgs.grid(column=0,row=0,padx=2,pady=2)

#this is the candidate image frame
f3 = LabelFrame(root, text="Similar Images",height=90*5,width=150*7)
f3.pack(fill="both", expand="yes")

root.mainloop()

print "cleaning up..."
for f in tempFiles:
    os.system("rm "+f)

print "goodbye"


