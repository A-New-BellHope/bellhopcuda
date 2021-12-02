#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: ./compare_shdfil.sh MunkB_Coh"
    echo "No path, no .shd"
    exit 1
fi

cxxfile="temp/cxx/$1.shd"
forfile="temp/FORTRAN/$1.shd"

if [ ! -f $cxxfile ] || [ ! -f $forfile ]; then
    echo "Could not find $cxxfile or $forfile"
    exit 1
fi

dd if=$cxxfile of=temp1.bin bs=$((2002*4)) count=10 2> /dev/null
dd if=$forfile of=temp2.bin bs=$((2002*4)) count=10 2> /dev/null

diff <(xxd temp1.bin) <(xxd temp2.bin)

rm -f temp1.bin temp2.bin
