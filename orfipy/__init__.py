#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 18:57:59 2020

@author: usingh
"""
import psutil
import multiprocessing

###Estimate procs and chunk
#total cores available
_max_procs=int(multiprocessing.cpu_count())
#total system memory in MB
mem_values=psutil.virtual_memory()
_max_mem=mem_values.total >> 20

"""
A long seq is ~ 150-200 MB (around chr1 on human)

roughly determine how many long seqs and their ORFs can fit in total_memory
"""
_long_seq_bytes=min(2e8,mem_values[0]/60)
#limit on how many long seqs to load at once; max of 3 gives best performance; more than 3 slows down
_long_limit=min(max(int(mem_values[0]/2e8/60),1),3)
