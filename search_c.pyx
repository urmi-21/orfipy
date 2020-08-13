#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:10:37 2020

@author: usingh
"""
import re
from libc.string cimport memcpy
from libc.stdlib cimport malloc, free



cdef int from_str_to_chararray(source, char *dest) except -1:
    cdef size_t source_len = len(source) 
    cdef bytes as_bytes = source.encode('ascii')    #hold reference to the underlying byte-object
    cdef char *as_ptr = <const char *>(as_bytes)
    memcpy(dest, as_ptr, source_len)
    return 0

def search(seq):
    #convert seq to array
    cdef int seq_len = len(seq) 
    cdef char *seq_array= <char *> malloc(seq_len*sizeof(char))
    from_str_to_chararray(seq,seq_array)
    
    #print("In array: ", seq_array)
    starts=['ATG']
    stops=['TAA','TAG','TGA']
    cdef int starts_found[3]
    starts_found[:]=[-1,-1,-1]#starts in frame 1 2 3
    cdef list orfs=[]
    cdef list allstarts=[]
    cdef int this_frame
    for i in range(seq_len-2):
        this_frame=i%3
        thiscodon=seq_array[i:i+3]
        #print(thiscodon)
        if thiscodon in starts:
            allstarts.append(i)
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
    print(len(allstarts))
    free(seq_array)
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
    print(len(start_positions))


