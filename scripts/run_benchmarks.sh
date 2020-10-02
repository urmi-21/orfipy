#!/bin/bash

#run benchmark
python benchmark.py $1

#test sequences match
python compare_fasta_files.py ../testdata/getorf_d ../testdata/orfm_d ../testdata/orfipy_testout/d
