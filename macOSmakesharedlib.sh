#!/bin/bash
# Author: Johannes Buchner (C) 2015
# For creating a shared library (libcuba.so).
# Edit for Mac OSX

sed 's/CFLAGS = -O3 -fomit-frame-pointer/CFLAGS = -O3 -fPIC -fomit-frame-pointer/g' --in-place makefile
echo "rebuilding libcuba.a archive"
make -B libcuba.a
echo "unpacking libcuba.a"
FILES=$(ar xv libcuba.a |sed 's/x - //g')
echo "making libcuba.so"
#echo gcc -shared -Wall $FILES -lm -o libcuba.so
#gcc -shared -Wall $FILES -lm -o libcuba.so
echo gcc -dynamiclib -Wall $FILES -lm -o libcuba.so
gcc -dynamiclib -Wall Vegas.o Vegas_.o llVegas.o llVegas_.o Suave.o Suave_.o llSuave.o llSuave_.o Divonne.o Divonne_.o llDivonne.o llDivonne_.o Cuhre.o Cuhre_.o llCuhre.o llCuhre_.o Fork.o Fork_.o Global.o Global_.o Data.o -lm -o libcuba.so

echo "created libcuba.so : now copy it into the location /usr/local/lib and pip install pymultinest for use with python"
rm $FILES

