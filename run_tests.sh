#!/bin/bash

#set -x

if [ -z $1 ]; then
    echo "Usage: ./run_tests.sh tl_match"
    exit 1
fi

ignore=0
if [[ $2 == "ignore" || $3 == "ignore" ]]; then
    ignore=1
fi

desiredresult=0
failedmsg="failed"
if [[ $2 == "shouldfail" ]]; then
    desiredresult=1
    failedmsg="failed to fail"
fi

check_fail () {
    i=$1
    if [[ $i != "0" ]]; then
        i=1
    fi
    if [[ $i != $desiredresult ]]; then
        echo "$3: $2 $failedmsg"
        if [[ $ignore != "1" ]]; then
            exit
        fi
    fi
}

tldesiredresult=0
tlfailedmsg="failed to match"
if [[ $2 == "tlshouldfail" ]]; then
    tldesiredresult=1
    tlfailedmsg="were not supposed to match, but they did"
fi

tl_check_fail () {
    i=$1
    if [[ $i != "0" ]]; then
        i=1
    fi
    if [[ $i != $tldesiredresult ]]; then
        echo "$2: TL results $tlfailedmsg"
        if [[ $ignore != "1" ]]; then
            exit
        fi
    fi
}

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
    runname="bellhopcxx single-threaded"
    echo $runname
    ./bin/bellhopcxx -1 test/cxx1/$1
    check_fail $? $runname $1
    runname="bellhopcxx multi-threaded"
    echo $runname
    ./bin/bellhopcxx test/cxxmulti/$1
    check_fail $? $runname $1
    runname="bellhopcuda"
    echo $runname
    ./bin/bellhopcuda test/cuda/$1
    check_fail $? $runname $1
    forres=0
    runname="BELLHOP"
    echo $runname
    forout=$(time ../bellhop/Bellhop/bellhop.exe test/FORTRAN/$1 2>&1)
    if [[ $? != "0" ]]; then forres=1; fi
    if [[ "$forout" == *"STOP Fatal Error"* ]]; then forres=1; fi
    check_fail $forres $runname
    python3 compare_shdfil.py $1
    tl_check_fail $? $1
}

while read -u 10 line || [[ -n $line ]]; do
    if [[ $line == //* ]]; then
        continue
    fi
    run_test $line
done 10<${1}.txt

echo "============================"
echo "Tests completed successfully"
echo "============================"
