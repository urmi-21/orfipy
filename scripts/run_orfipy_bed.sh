#!/bin/bash

for i in $(seq $4)
do
    orfipy --min $3 --bed orfipy_bed --outdir $2 --between-stop $1
done
