#!/bin/bash

for i in $(seq $4)
do
    orfipy --min $3 --dna orfipy_d --pep orfipy_p --outdir $2 --between-stop $1
done
