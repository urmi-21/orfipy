#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:43:24 2020

@author: usingh

Use this script to benchmark runtimes of different tools

arguments[1]<-input_fasta_file
"""

from pyrpipe import pyrpipe_engine as pe
import sys

N="5"
minlen="3"
outdir="times_out"
testseq=sys.argv[1]
#compare orfipy, orfm, getorf
orfipy_cmd=["bash", "run_orfipy.sh",testseq,outdir,minlen,N]

orfm_cmd=["bash", "run_orfm.sh",testseq,outdir,minlen,N]

getorf_cmd=["bash", "run_getorf.sh",testseq,outdir,minlen,N]

#for i in range(N):
pe.execute_command(getorf_cmd,objectid=testseq,command_name="getorf")
pe.execute_command(orfipy_cmd,objectid=testseq,command_name="orfipy")
pe.execute_command(orfm_cmd,objectid=testseq,command_name="orfm")



#compare seqs
#python compare_fasta_files.py ../testdata/getorf_d ../testdata/orfm_d ../testdata/orfipy_testout/d
#python compare_fasta_files.py ../testdata/getorf_p ../testdata/orfm_p ../testdata/orfipy_testout/p


###Compare orfipy and getorf -3 option
#orfipy --min 3 --dna 3_d --pep 3_p --outdir orfipy_testout --start ATG --partial-3 testseq.fa
#getorf -find 3 -min 3 -outseq getorf3_d -sequence testseq.fa
#getorf -find 1 -min 3 -outseq getorf3_p -sequence testseq.fa

#python compare_fasta_files.py ../testdata/getorf3_p ../testdata/orfipy_testout/3_p
#python compare_fasta_files.py ../testdata/getorf3_d ../testdata/orfipy_testout/3_d



