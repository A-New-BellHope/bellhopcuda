#!/bin/bash

#set -x

run_test () {
    echo $1
    rm -f temp/cxx/$1.*
    rm -f temp/FORTRAN/$1.*
    cp temp/in/$1.* temp/cxx/
    cp temp/in/$1.* temp/FORTRAN/
    echo "bellhopcxx"
    ./bin/bellhopcxx -1 temp/cxx/$1
    #./bin/bellhopcxx temp/cxx/$1
    cxxres=$?
    forres=0
    echo "BELLHOP"
    forout=$(../thirdparty/at_2020_11_4/Bellhop/bellhop.exe temp/FORTRAN/$1 2>&1)
    if [[ $? != "0" ]]; then forres=1; fi
    if [[ "$forout" == *"STOP Fatal Error"* ]]; then forres=1; fi
}

while read -u 10 line || [[ -n $line ]]; do
    if [[ $line == //* ]]; then
        continue
    fi
    run_test $line
    if [[ $cxxres != "0" ]]; then
        echo "$line: bellhopcxx failed"
    fi
    if [[ $forres != "0" ]]; then
        echo "$line: BELLHOP failed"
    fi
    if [[ $cxxres != "0" || $forres != "0" ]]; then
        exit 1
    fi
done 10<tests_pass.txt

while read -u 11 line || [[ -n $line ]]; do
    if [[ $line == //* ]]; then
        continue
    fi
    run_test $line
    if [[ $cxxres == "0" ]]; then
        echo "$line: bellhopcxx failed to fail"
    fi
    if [[ $forres == "0" ]]; then
        echo "$line: BELLHOP failed to fail"
    fi
    if [[ $cxxres == "0" || $forres == "0" ]]; then
        exit 1
    fi
done 11<tests_fail.txt

echo "Tests completed"
