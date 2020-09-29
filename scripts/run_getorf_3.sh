#!/bin/bash
for i in $(seq $4)
do
    getorf -find 1 -min $3 -outseq "$2/getorf_3_p" -sequence $1
    getorf -find 3 -min $3 -outseq "$2/getorf_3_d" -sequence $1
done
