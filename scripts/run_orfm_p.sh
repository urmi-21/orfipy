#!/bin/bash
for i in $(seq $4)
do
    orfm -m $3 $1 > "$2/orfm_p_only"
done
