#!/bin/bash

#set -x

if [ -z $1 ]; then
    echo "Usage: ./run_tests.sh ray_tests"
    exit 1
fi

run_test () {
    echo $1
    if [ ! -f "test/in/$1.env" ]; then
        echo "test/in/$1.env does not exist"
        exit 1
    fi
    mkdir -p test/cxx1
    mkdir -p test/cxxmulti
    mkdir -p test/cuda
    mkdir -p test/FORTRAN
    rm -f test/cxx1/$1.*
    rm -f test/cxxmulti/$1.*
    rm -f test/cuda/$1.*
    rm -f test/FORTRAN/$1.*
    cp test/in/$1.* test/cxx1/
    cp test/in/$1.* test/cxxmulti/
    cp test/in/$1.* test/cuda/
    cp test/in/$1.* test/FORTRAN/
    echo "bellhopcxx single-threaded"
    ./bin/bellhopcxx -1 test/cxx1/$1
    cxx1res=$?
    echo "bellhopcxx multi-threaded"
    ./bin/bellhopcxx test/cxxmulti/$1
    cxxmultires=$?
    echo "bellhopcuda"
    ./bin/bellhopcuda test/cuda/$1
    cudares=$?
    forres=0
    echo "BELLHOP"
    forout=$(time ../bellhop/Bellhop/bellhop.exe test/FORTRAN/$1 2>&1)
    if [[ $? != "0" ]]; then forres=1; fi
    if [[ "$forout" == *"STOP Fatal Error"* ]]; then forres=1; fi
}

while read -u 10 line || [[ -n $line ]]; do
    if [[ $line == //* ]]; then
        continue
    fi
    run_test $line
    if [[ $cxx1res != "0" ]]; then
        echo "$line: bellhopcxx single-threaded failed"
    fi
    if [[ $cxxmultires != "0" ]]; then
        echo "$line: bellhopcxx multi-threaded failed"
    fi
    if [[ $cudares != "0" ]]; then
        echo "$line: bellhopcuda failed"
    fi
    if [[ $forres != "0" ]]; then
        echo "$line: BELLHOP failed"
    fi
    if [[ $cxx1res != "0" || $cxxmultires != "0" || $cudares != "0" || $forres != "0" ]]; then
        exit 1
    fi
done 10<${1}_pass.txt

while read -u 11 line || [[ -n $line ]]; do
    if [[ $line == //* ]]; then
        continue
    fi
    run_test $line
    if [[ $cxx1res == "0" ]]; then
        echo "$line: bellhopcxx single-threaded failed to fail"
    fi
    if [[ $cxxmultires == "0" ]]; then
        echo "$line: bellhopcxx multi-threaded failed to fail"
    fi
    if [[ $cudares == "0" ]]; then
        echo "$line: bellhopcuda failed to fail"
    fi
    if [[ $forres == "0" ]]; then
        echo "$line: BELLHOP failed to fail"
    fi
    if [[ $cxx1res == "0" || $cxxmultires == "0" || $cudares == "0" || $forres == "0" ]]; then
        exit 1
    fi
done 11<${1}_fail.txt

echo "============================"
echo "Tests completed successfully"
echo "============================"
