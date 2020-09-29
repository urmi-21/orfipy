#!/bin/bash
for i in $(seq $4)
do
    orfm -m $3 -t "$2/orfm_d" $1 > "$2/orfm_p"
done
