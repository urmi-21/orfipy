#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:43:24 2020

@author: usingh

Use this script to benchmark runtimes of different tools

args[1]<-input_fasta_file
args[2]<-out dir
[3]<-Num trials
[4]<-Min len of ORF should be multiple of 3 for orfm to work
"""

from pyrpipe import pyrpipe_engine as pe
import sys
import os

#print(sys.argv)
outdir=sys.argv[2]
if not os.path.exists(outdir):
    os.makedirs(outdir)

N=sys.argv[3]
minlen=sys.argv[4]

testseq=sys.argv[1]
#compare orfipy, orfm, getorf
orfipy_cmd_b=["bash", "run_orfipy_bed.sh",testseq,outdir,minlen,N]
orfipy_cmd=["bash", "run_orfipy.sh",testseq,outdir,minlen,N]
orfm_cmd=["bash", "run_orfm.sh",testseq,outdir,minlen,N]
orfm_cmd_p=["bash", "run_orfm_p.sh",testseq,outdir,minlen,N]
getorf_cmd=["bash", "run_getorf.sh",testseq,outdir,minlen,N]

pe.execute_command(orfipy_cmd_b,objectid="test",command_name="orfipy_b")
pe.execute_command(orfipy_cmd,objectid="test",command_name="orfipy")
pe.execute_command(orfm_cmd_p,objectid="test",command_name="orfm_p")
pe.execute_command(orfm_cmd,objectid="test",command_name="orfm")
pe.execute_command(getorf_cmd,objectid="test",command_name="getorf")


###Compare orfipy and getorf -3 option
#orfipy_cmd=["bash", "run_orfipy_3.sh",testseq,outdir,minlen,N]
#getorf_cmd=["bash", "run_getorf_3.sh",testseq,outdir,minlen,N]
#pe.execute_command(getorf_cmd,objectid="test",command_name="getorf_3")
#pe.execute_command(orfipy_cmd,objectid="test",command_name="orfipy_3")

