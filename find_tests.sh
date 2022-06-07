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

testsdir="../bellhop/tests"

echo "3D ray:"
find $testsdir -name "*.env" | while read f; do
    res=$(sed -n '10,$s|^'\''R[^/]*3'\''|&|p' $f)
    if [ -n "$res" ]; then
        echo "$f: $res"
    fi
done

echo "Nx2D ray:"
find $testsdir -name "*.env" | while read f; do
    res=$(sed -n '10,$s|^'\''R[^/]*2'\''|&|p' $f)
    if [ -n "$res" ]; then
        echo "$f: $res"
    fi
done

exit

echo "Arrivals:"
find $testsdir -name "*.env" | while read f; do
    res=$(sed -n '10,$s|^'\''[Aa][^3/~]*'\''[^0]*$|&|p' $f)
    if [ -n "$res" ]; then
        echo "$f: $res"
    fi
done

exit

echo "Eigenray:"
find $testsdir -name "*.env" | while read f; do
    res=$(sed -n '10,$s|^'\''E[^3/]*'\''|&|p' $f)
    if [ -n "$res" ]; then
        echo "$f: $res"
    fi
done

exit

echo "Ray:"
find $testsdir -name "*.env" | while read f; do
    res=$(sed -n '10,$s|^'\''R[^3/]*'\''|&|p' $f)
    if [ -n "$res" ]; then
        echo "$f: $res"
    fi
done

echo "Coherent:"
find $testsdir -name "*.env" | while read f; do
    res=$(sed -n '10,$s|^'\''C[CRSbBgG][^3/]*'\''|&|p' $f)
    if [ -n "$res" ]; then
        echo "$f: $res"
    fi
done

echo "Semi-coherent:"
find $testsdir -name "*.env" | while read f; do
    res=$(sed -n '10,$s|^'\''S[CRSbBgG][^3/]*'\''|&|p' $f)
    if [ -n "$res" ]; then
        echo "$f: $res"
    fi
done

echo "Incoherent:"
find $testsdir -name "*.env" | while read f; do
    res=$(sed -n '10,$s|^'\''I[CRSbBgG][^3/]*'\''|&|p' $f)
    if [ -n "$res" ]; then
        echo "$f: $res"
    fi
done
