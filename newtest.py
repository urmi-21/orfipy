#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:07:13 2020

@author: usingh
"""
import sys
import time
import re
import pyximport; pyximport.install()
import search_c
import orfipy_core as oc

def search(seq):
    l=len(seq)
    starts=['ATG']
    stops=['TAA','TAG','TGA']
    starts_found=[-1,-1,-1]#starts in frame 1 2 3
    orfs=[]
    for i in range(l-2):
        this_frame=i%3
        thiscodon=seq[i:i+3]
        #print(thiscodon)
        if thiscodon in starts:
            #print('strt',i)
            #add only if no prev start, if a start is already in list then ignore
            if starts_found[this_frame]<0:
                starts_found[this_frame]=i
        elif thiscodon in stops:
            #print('stp',i)
            #print(starts_found)
            #if a start exist
            if starts_found[this_frame]>=0:
                orfs.append((starts_found[this_frame],i+2))
                #print(seq[starts_found[this_frame]:i+2])
                starts_found[this_frame]=-1
    return orfs

def research(seq):
    starts=['ATG']
    stops=['TAA','TAG','TGA']
    start_positions=[]
    stop_positions=[]
    for c in starts:
        start_positions.extend([m.start() for m in re.finditer(c,seq)])
    for c in stops:
        stop_positions.extend([m.start() for m in re.finditer(c,seq)])
    #print()
    
alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
def reverse_complement(seq):    
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

test='ATGAATGAACCCTAGATAGCC'
with open(sys.argv[1]) as f:
    d=f.read().splitlines()

d.pop(0)#remove header
test=''.join(d)
print(len(test))

start_time = time.time()
duration=time.time() - start_time
oc.get_orfs_new(test,'test',30)
print('XXXXXXXXXXXXXXXXXXXXXXXX')

oc.get_orfs_new(reverse_complement(test),'testrc',30,rev_com=True)
print('t1',duration)


sys.exit(0)
print('XXXXXXXXXXXXXXXXXXXXXXXX')
start_time = time.time()
#search(test)
duration=time.time() - start_time
print('t1',duration)

start_time=time.time()
research(test)
duration=time.time() - start_time
print('t2',duration)

start_time=time.time()
search_c.search(test)
duration=time.time() - start_time
print('t3',duration)

start_time=time.time()
search_c.research(test)
duration=time.time() - start_time
print('t3',duration)