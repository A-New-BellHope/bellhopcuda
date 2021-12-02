#!/bin/bash

find ../thirdparty/at_2020_11_4/tests -name "*.env" | while read f; do
    res=$(sed -n '10,$s/'\''C[CRSbBgG][^3]*'\''/&/p' $f)
    if [ -n "$res" ]; then
        echo "$f: $res"
    fi
done
