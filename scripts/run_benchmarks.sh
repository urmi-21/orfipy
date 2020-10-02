#!/bin/bash
set -e
#infiles 
#$1 <- fasta file
#$2 <- outdir
#$3 <- num trials
#$4<- minlen orf orfs (multiple of 3)
#example bash run_benchmarks.sh ../testdata/testseq.fa b1 5 30

#run benchmark
python benchmark.py $1 $2 $3 $4

#test all sequences match
echo "Comparing nucleotide sequences"
python compare_fasta_files.py $2/getorf_d $2/orfm_d $2/orfipy_d
echo "Comparing peptide sequences"
python compare_fasta_files.py $2/getorf_p $2/orfm_p $2/orfipy_p

#get pyrpipelog file
plog=$(ls -Art pyrpipe_logs/*pyrpipe.log | tail -n 1)

#summarize runtimes
pyrpipe_diagnostic.py benchmark -t $2/tmp $plog

timecsv="$2/tmp/benchmark_reports/time_per_program.csv"

#print avg time for each program
echo " "

head -1 $timecsv | awk -F "," '{print $1 "\tAvg time(sec)"}' |tr "," "\t"
sed '1d' $timecsv | awk -v N=$3 -F "," '{printf "%s,%4.2f\n", $1 , $2/N}' | tr "," "\t"

