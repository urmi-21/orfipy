#!/bin/bash

#example: bash basic_benchmark_fq.sh ../testdata/bmdata/SRR976159_1.fastq.gz 66 fqgz_out

infile=$1
min=$2
outdir=$3

mkdir -p $3
time orfipy --min $2 --bed b --outdir $3 --between-stops $1
time orfipy --min $2 --dna d --pep p --outdir $3 --between-stops $1
time orfm -m $2 -t "$3/orfm_d" $1 > "$3/orfm_p"
time orfm -m $2 $1 > "$3/orfm_p_only"

