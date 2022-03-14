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

if [ -z "$1" ]; then
    echo "Usage: ./compare_shdfil.sh MunkB_Coh"
    echo "No path, no .shd"
    exit 1
fi

cxxfile="temp/cxx/$1.shd"
forfile="temp/FORTRAN/$1.shd"
xxdopts="-e"

compare () {
    xxdopts="-e -c 8"
    mkfifo cxxpipe
    mkfifo forpipe
    xxd $xxdopts $1 >cxxpipe &
    xxd $xxdopts $2 >forpipe &
    
    while true; do
        if ! read -u 10 cxxline && [[ -z $cxxline ]]; then
            break
        fi
        if ! read -u 11 forline && [[ -z $forline ]]; then
            break
        fi
        cxxaddr=$(echo "$cxxline" | sed -r "s/(\w+): .*/\1/")
        foraddr=$(echo "$forline" | sed -r "s/(\w+): .*/\1/")
        if [[ "$cxxaddr" != "$foraddr" ]]; then
            echo "Addresses did not match: $cxxaddr $foraddr"
            exit 1
        fi
        cxxdata=$(echo "$cxxline" | sed -r "s/\w+: (\w+ \w+) .*/\1/")
        fordata=$(echo "$forline" | sed -r "s/\w+: (\w+ \w+) .*/\1/")
        if [[ "$cxxdata" != "$fordata" ]]; then
            echo "$cxxaddr: CXX $cxxdata  FOR $fordata"
        fi    
    done 10<cxxpipe 11<forpipe
    
    rm cxxpipe
    rm forpipe
}

if [[ $2 == "full" ]]; then
    #compare $cxxfile $forfile
    diff <(xxd $xxdopts $cxxfile) <(xxd $xxdopts $forfile)
    exit 0
elif [[ -n "$2" ]]; then
    echo "Invalid second argument, must be \"full\" or omitted"
    exit 1
fi

if [ ! -f $cxxfile ] || [ ! -f $forfile ]; then
    echo "Could not find $cxxfile or $forfile"
    exit 1
fi

dd if=$forfile of=temp1.bin bs=1 count=4 2> /dev/null
reclen0=$(xxd temp1.bin)
rm -f temp1.bin
reclen1=$(echo "$reclen0" | sed -r "s/00000000: (\w\w)(\w\w) (\w\w)(\w\w).*/\4\3\2\1/")
echo "$reclen1"
reclen2=$((0x$reclen1*4))

dd if=$cxxfile of=temp1.bin bs=$reclen2 count=10 2> /dev/null
dd if=$forfile of=temp2.bin bs=$reclen2 count=10 2> /dev/null

#compare temp1.bin temp2.bin
diff <(xxd $xxdopts temp1.bin) <(xxd $xxdopts temp2.bin)

rm -f temp1.bin temp2.bin
