#!/bin/bash

#set -x

if [[ -z $1 || -z $2 ]]; then
    echo "Usage: ./run_tests.sh (ray/tl/eigen/arr) tests_list [(nothing)/shouldfail/shouldnotmatch] [ignore]"
    exit 1
fi

if [[ $1 != "ray" && $1 != "tl" && $1 != "eigen" && $1 != "arr" ]]; then
    echo "Invalid run type (first arg), must be: ray, tl, eigen, arr"
    exit 1
fi
runtype=$1

ignore=0
if [[ $3 == "ignore" || $4 == "ignore" ]]; then
    ignore=1
fi

desiredresult=0
failedmsg="failed"
if [[ $3 == "shouldfail" ]]; then
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

mdesiredresult=0
mfailedmsg="failed to match"
if [[ $3 == "shouldnotmatch" ]]; then
    mdesiredresult=1
    mfailedmsg="were not supposed to match, but they did"
fi

m_check_fail () {
    i=$1
    if [[ $i != "0" ]]; then
        i=1
    fi
    if [[ $i != $mdesiredresult ]]; then
        echo "$2: $runtype results $mfailedmsg"
        if [[ $ignore != "1" ]]; then
            exit
        fi
    fi
}

run_test () {
    echo ""
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
    forres=0
    runname="BELLHOP"
    echo $runname
    forout=$(time ../bellhop/Bellhop/bellhop.exe test/FORTRAN/$1 2>&1)
    if [[ $? != "0" ]]; then forres=1; fi
    if [[ "$forout" == *"STOP Fatal Error"* ]]; then forres=1; fi
    check_fail $forres $runname
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
    if [[ $runtype == "tl" ]]; then
        python3 compare_shdfil.py $1
        m_check_fail $? $1
    elif [[ $runtype == "arr" ]]; then
        python3 compare_arrivals.py $1
        m_check_fail $? $1
    else
        echo "Automated checker not installed for $runtype runs"
    fi
}

while read -u 10 line || [[ -n $line ]]; do
    if [[ $line == //* ]]; then
        continue
    fi
    run_test $line
done 10<$2.txt

echo "============================"
echo "Tests completed successfully"
echo "============================"
