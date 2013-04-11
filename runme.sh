#!/bin/bash
#convert to pgms
echo "converting jpegs to pgms and storing in pgm folder"
cd pictures
for i in `ls *.jpg`; do convert $i $i.pgm; done
mkdir ../pgm
mv *.pgm ../pgm
cd ..
cp sift pgm/
cd pgm
echo "computing sift vectors for pgms"
#compute sift vectors
for i in `ls *.pgm`;
do
    if ls $i.key
    then
        echo "here"
    else
        ./sift <$i> $i.key;
    fi
done

cd ..
echo "creating sift leech lattice hashed database"
#make the database
python knn2.py
echo "dirtying up the images and converting to pgm"
#dirty up the images with some blur and compression
cd pictures
for i in `ls *.jpg`; 
do NUMBER=$[ ((( $RANDOM % 4 )  + 1 ) *90)]
convert -quality 20 -rotate $NUMBER -adaptive-blur 2x2 $i $i.pgm  
done
mkdir ../pgmdirty
mv *.pgm ../pgmdirty
cd ..
cp sift pgmdirty/
cd pgmdirty
echo "matching dirtied images to original images"
#compute sift vectors
for i in `ls *.pgm`;
do
    if ls $i.key
    then
        echo "here"
    else
        ./sift <$i> $i.key;
    fi
done

cd ..
python test.py > output


