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
#print('XXXX',_max_mem,_max_procs)

"""
A long seq is ~ 4e8 bytes (around chr1 on human)

roughly determine how many long seqs and their ORFs can fit in total_memory/2
"""
_long_seq_bytes=4e8
_long_limit=max(min(int(mem_values[0]/_long_seq_bytes/30),5),1)

print('LONG limit',_long_limit)