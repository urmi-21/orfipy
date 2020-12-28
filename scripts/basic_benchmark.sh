#!/bin/bash

#example: bash basic_benchmark.sh ../testdata/bmdata/human.fa 99 human_basic_out

infile=$1
min=$2
outdir=$3

mkdir -p $3
time orfipy --min $2 --bed b --outdir $3 --between-stops $1
time orfipy --min $2 --dna d --pep p --outdir $3 --between-stops $1
time orfm -m $2 -t "$3/orfm_d" $1 > "$3/orfm_p"
time orfm -m $2 $1 > "$3/orfm_p_only"
time ( getorf -find 0 -min $2 -outseq "$3/getorf_p" -sequence $1 ; getorf -find 2 -min $2 -outseq "$3/getorf_d" -sequence $1 )

