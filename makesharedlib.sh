#!/bin/bash
# Author: Johannes Buchner (C) 2015
# For creating a shared library (libcuba.so).

sed 's/CFLAGS = -O3 -fomit-frame-pointer/CFLAGS = -O3 -fPIC -fomit-frame-pointer/g' --in-place makefile
echo "rebuilding libcuba.a archive"
make -B libcuba.a
echo "unpacking libcuba.a"
FILES=$(ar xv libcuba.a |sed 's/x - //g'|grep -v ' __\.SYMDEF SORTED__')
echo "making libcuba.so"
echo gcc -shared -Wall -lm $FILES -o libcuba.so
gcc -shared -Wall -lm $FILES -o libcuba.so
rm $FILES


