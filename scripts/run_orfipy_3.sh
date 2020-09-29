#!/bin/bash

for i in $(seq $4)
do
    orfipy --min $3 --dna orfipy_3_d --pep orfipy_3_p --outdir $2 --start ATG --partial-3 $1
done
