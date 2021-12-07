#!/bin/bash

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
