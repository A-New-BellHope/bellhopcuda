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
    echo "Usage: ./run_tests.sh (ray/tl/eigen/arr)(3d) tests_list [(nothing)/shouldfail/shouldnotmatch] [ignore]"
    exit 1
fi

runtype=`echo $1 | sed 's/3d//g'`
threedopt="-2"
bhexec=bellhop.exe
if [[ $1 == *3d ]]; then
    threedopt="-3"
    bhexec=bellhop3d.exe
fi
if [[ $runtype == "ray" || $runtype == "eigen" ]]; then
    comparepy=compare_ray.py
elif [[ $runtype == "tl" ]]; then
    comparepy=compare_shdfil.py
elif [[ $runtype == "arr" ]]; then
    comparepy=compare_arrivals.py
else
    echo "$runtype is not a valid run type (valid types are ray/tl/eigen/arr, optionally + 3d)"
fi

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
finalmsg="passed"
if [[ $3 == "shouldfail" ]]; then
    desiredresult=1
    failedmsg="failed to fail"
    finalmsg="failed"
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
        echo "$3: $2 results $mfailedmsg"
        if [[ $ignore != "1" ]]; then
            return 1
        fi
    fi
    return 0
}

dotexe=""
if [ -f ./bin/bellhopcxx.exe ]; then
    dotexe=".exe"
fi

compare_results () {
    runname="$1"
    dir=$2
    envfil=$3
    if [[ $desiredresult == "1" ]]; then
        echo "Skipping results comparison because in shouldfail mode"
        return 0
    elif [[ $runtype != "tl" && ( $dir == "cxxmulti" || $dir == "cuda" ) ]]; then
        echo "Skipping $runtype results comparison for $dir"
        return 0
    fi
    python3 $comparepy $envfil $dir
    m_check_fail $? "$runname" $envfil
}

run_fortran () {
    envfil=$1
    runname="BELLHOP"
    echo $runname
    forres=0
    forout=$(time ../bellhop/Bellhop/$bhexec test/FORTRAN/$envfil 2>&1)
    if [[ $? != "0" ]]; then forres=1; fi
    if [[ "$forout" == *"STOP Fatal Error"* ]]; then forres=1; fi
    check_fail $forres "$runname" $envfil
}

run_cxx1 () {
    envfil=$1
    runname="bellhopcxx single-threaded"
    echo $runname
    ./bin/bellhopcxx$dotexe -1 $threedopt test/cxx1/$envfil
    check_fail $? "$runname" $envfil
}

run_check_and_compare () {
    runname="$1"
    prog=$2
    dir=$3
    envfil=$4
    echo $runname
    $prog $threedopt test/$dir/$envfil
    check_fail $? "$runname" $envfil || return 1
    compare_results "$runname" $dir $envfil || return 2
    return 0
}

try_a_few_times () {
    runname="$1"
    prog=$2
    dir=$3
    envfil=$4
    count=1
    maxcount=5
    while true; do
        run_check_and_compare "$runname" $prog $dir $envfil
        res=$?
        if [[ $res < 2 ]]; then return $res; fi
        count=$((count + 1))
        if [ "$count" -gt "$maxcount" ]; then
            echo "Failed to match after $maxcount re-rolls"
            return 1
        fi
        echo "Re-rolling RNG $count of $maxcount..."
    done
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
    run_fortran $1 &
    fortranpid=$!
    run_cxx1 $1 &
    cxx1pid=$!
    wait $fortranpid || exit 1
    wait $cxx1pid || exit 1
    compare_results "bellhopcxx single-threaded" cxx1 $1 || exit 1
    try_a_few_times "bellhopcxx multi-threaded" ./bin/bellhopcxx$dotexe cxxmulti $1 || exit 1
    if [ -f ./bin/bellhopcuda$dotexe ]; then
        try_a_few_times "bellhopcuda" ./bin/bellhopcuda$dotexe cuda $1 || exit 1
    else
        echo "bellhopcuda$dotexe not found ... ignoring"
    fi
}

while read -u 10 line || [[ -n $line ]]; do
    if [[ $line == //* ]]; then
        continue
    fi
    run_test $line
done 10<$2.txt

echo ""
echo "============================="
echo "All tests $finalmsg successfully"
echo "============================="
