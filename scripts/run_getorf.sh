#!/bin/bash
for i in $(seq $4)
do
    getorf -find 0 -min $3 -outseq "$2/getorf_p" -sequence $1
    getorf -find 2 -min $3 -outseq "$2/getorf_d" -sequence $1
done
