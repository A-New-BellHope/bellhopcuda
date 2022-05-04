#!/bin/bash
# bellhopcxx / bellhopcuda - C++/CUDA port of BELLHOP underwater acoustics simulator
# Copyright (C) 2021-2022 The Regents of the University of California
# c/o Jules Jaffe team at SIO / UCSD, jjaffe@ucsd.edu
# Based on BELLHOP, which is Copyright (C) 1983-2020 Michael B. Porter
# 
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

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

if [[ $2 == *.txt ]]; then
    echo "BELLHOP syntax prohibits including the file extension on input files (drop the .txt)"
    exit 1
fi

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
    if [ -f ./bin/bellhopcuda ]; then
        ./bin/bellhopcuda test/cuda/$1
        check_fail $? $runname $1
    else
        echo "bellhopcuda not found ... ignoring"
    fi
    if [[ $runtype == "ray" ]]; then
        python3 compare_ray.py $1
        m_check_fail $? $1
    elif [[ $runtype == "tl" ]]; then
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
