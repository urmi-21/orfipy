#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 20:01:56 2020

@author: usingh
"""
import os
import sys
import time
import multiprocessing
from contextlib import closing
import gc
import psutil
import pyximport; pyximport.install()
import orfipy.orfipy_core as oc
from pyfaidx import Fasta





def worker(thisseq,thisseq_rc,thisname,minlen,strand,starts,stops,bed12,bed,dna,rna,pep):
    """
    start worker
    """
    #call orf function
    res=oc.find_orfs(thisseq,thisseq_rc,thisname,minlen,strand,starts,stops,bed12,bed,dna,rna,pep)
    return res

def worker2(arglist):
    """
    start worker
    """
    #call orf function
    res=oc.find_orfs(arglist[0],arglist[1],arglist[2],arglist[3],arglist[4],arglist[5],arglist[6],arglist[7],arglist[8],arglist[9],arglist[10],arglist[11])
    return res


def start_multiprocessor(seqs,minlen,procs,chunk_size,strand,starts,stops,bed12,bed,dna,rna,pep):
    #process = psutil.Process(os.getpid())
    #poolargs contain data to be passed to mp workers
    poolargs=[]
    #Use a memory limit to roughly limit memory usage to  order of _MEMLIMIT
    #this is helpful in making code run on low mwmory systems and with python 3.8 and lower if data is more than 2GB in size
    _MEMLIMIT=chunk_size*1000000 #convert MB to bytes
    #total bytes read in memory
    total_read_bytes=0
    #all read bytes
    cummulative_read_bytes=0
    
    #process sequences in fasta file
    for s in list(seqs.keys()):
        
        thisname=s
        thisseq=str(seqs[s])
        #ignore if seq is < minlen
        if len(thisseq)<minlen:
            continue
        #add bytes read
        this_read=len(thisseq.encode('utf-8'))
        thisseq_rc=None
        if strand == 'b' or strand =='r':
            thisseq_rc=seqs[s][:].complement.reverse.seq
            #add bytes read of rev_com seq
            this_read=this_read*2
        
        #add read bytes to total_read_bytes
        total_read_bytes+=this_read
        cummulative_read_bytes+=this_read
        #print(s,total_read_bytes,this_read)
        
        #add to poolargs; if limit is reached this will be reset
        poolargs.append([thisseq,thisseq_rc,thisname,minlen,strand,starts,stops,bed12,bed,dna,rna,pep])
        
        #if total_read_bytes is more than memory limit
        if total_read_bytes+1000000 >= _MEMLIMIT:
            #mem=float(process.memory_info().rss)/1000000.0
            #print('Processed {0:.2f} MB:'.format(cummulative_read_bytes/1000000),'this c',total_read_bytes,'Memory usage: {0:.2f}MB'.format(mem), end="\r", flush=True,file=sys.stderr )
            #print('Processing {0:d} bytes:'.format(cummulative_read_bytes),'Current Memory usage: {0:.2f}MB'.format(mem), end="\r", flush=True,file=sys.stderr )
            print('Processed {0:d} bytes'.format(cummulative_read_bytes), end="\r", flush=True,file=sys.stderr)
            #process seqs that were added to poolargs
            with closing( multiprocessing.Pool(procs) ) as p:
                results_inner = p.imap_unordered(worker2, poolargs, 100)
            #print('#####')
            #convert results_inner from generator object to a list
            rlist=list(results_inner)
            #write these results to file
            write_results(rlist,bed12,bed,dna,rna,pep)
            #perfor GC
            del results_inner
            del poolargs
            gc.collect()
            #create empty list
            poolargs=[]
            #reset total read bytes for next batch
            total_read_bytes=0
            
            
        
        #after loop poolargs contains seq; process these 
    #print('executing outer',len(poolargs))
    
    if len(poolargs) > 0:
        print('Processed {0:d} bytes'.format(cummulative_read_bytes), end="\r", flush=True,file=sys.stderr)
        with closing( multiprocessing.Pool(procs) ) as p:
            results_inner = p.imap_unordered(worker2, poolargs, 100)
        #convert results_inner from generator object to a list
        rlist=list(results_inner)
        #print('###Got list##')
        #print(rlist)
        #write results to file
        write_results(rlist,bed12,bed,dna,rna,pep)
        #perform GC
        #print('start GC')
        del results_inner
        del poolargs
        gc.collect()
    
    print()
    

def worker_single(seqs,minlen,procs,strand,starts,stops,bed12,bed,dna,rna,pep):
    """
    Perform sequential processing

    Parameters
    ----------
    infasta : TYPE
        DESCRIPTION.
    minlen : TYPE
        DESCRIPTION.
    procs : TYPE
        DESCRIPTION.
    strand : TYPE
        DESCRIPTION.
    starts : TYPE
        DESCRIPTION.
    stops : TYPE
        DESCRIPTION.
    bed12 : TYPE
        DESCRIPTION.
    bed : TYPE
        DESCRIPTION.
    dna : TYPE
        DESCRIPTION.
    rna : TYPE
        DESCRIPTION.
    pep : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    print('orfipy single-mode',file=sys.stderr)
    results=[]
    for s in list(seqs.keys()):
        thisname=s
        thisseq=str(seqs[s])
        #ignore if seq is < minlen
        if len(thisseq)<minlen:
            continue
        thisseq_rc=None
        if strand == 'b' or strand =='r':
            thisseq_rc=seqs[s][:].complement.reverse.seq
        res=oc.find_orfs(thisseq,thisseq_rc,thisname,minlen,strand,starts,stops,bed12,bed,dna,rna,pep)
        results.append(res)
    
    write_results(results, bed12, bed, dna, rna, pep)
    




def write_results(results,bed12,bed,dna,rna,pep):
        
    #if no out type is specified
    if not (bed12 or bed or dna or rna or pep):
        #print('stdout')
        for reslist in results:
            print(reslist[0]) #reslist[0] contains bed
    else:
        if bed:
            bf=open(bed,'a')
        if bed12:
            b12f=open(bed12,'a')
        if dna:
            dnaf=open(dna,'a')
        if rna:
            rnaf=open(rna,'a')
        if pep:
            pepf=open(pep,'a')
        for reslist in results:
            if bed:
                bf.write(reslist[0]+'\n')
            if bed12:
                b12f.write(reslist[1]+'\n')
            if dna:
                dnaf.write(reslist[2]+'\n')
            if rna:
                rnaf.write(reslist[3]+'\n')
            if pep:
                pepf.write(reslist[4]+'\n')
        #close files
        if bed:
            bf.close()
        if bed12:
            b12f.close()
        if dna:
            dnaf.close()
        if rna:
            rnaf.close()
        if pep:
            pepf.close()
    
    
def init_result_files(bed12,bed,dna,rna,pep):
    #create empty files to append later
    if bed:
        f=open(bed,'w')
        f.close()
    if bed12:
        f=open(bed12,'w')
        f.close()
    if dna:
        f=open(dna,'w')
        f.close()
    if rna:
        f=open(rna,'w')
        f.close()
    if pep:
        f=open(pep,'w')
        f.close()

##########main################
def main(infasta,minlen,procs,single_mode,chunk_size,strand,starts,stops,bed12,bed,dna,rna,pep):
    init_result_files(bed12, bed, dna, rna, pep)    
    print("orfipy")
    ##start time
    start = time.time()
    
    if not procs:
        procs=multiprocessing.cpu_count()
    
    #estimate chunk_size
    if not chunk_size:
        #total mem in MB
        total_mem_MB=psutil.virtual_memory()[0]/1000000
        chunk_size=int(total_mem_MB/2)
    else:
        chunk_size=int(chunk_size)
    #check py < 3.8; if yes max chunk size can be 2000 other wise error is reported
    if sys.version_info[2] < 8 and chunk_size > 2000:
        chunk_size = 1900
                
    print("Setting chunk size {} MB".format(chunk_size),file=sys.stderr)
    #sys.exit(1)
    #read fasta file
    seqs = Fasta(infasta)
    #totalseqs=len(seqs.keys())
    
    
    if single_mode:
        worker_single(seqs, minlen, procs, strand, starts, stops, bed12, bed, dna, rna, pep)
        duration = time.time() - start
        print("Processed {0:d} sequences in {1:.2f} seconds".format(len(seqs.keys()),duration),file=sys.stderr)
        return
    else:
        start_multiprocessor(seqs, minlen, procs, chunk_size, strand, starts, stops, bed12, bed, dna, rna, pep)
        duration = time.time() - start
        print("Processed {0:d} sequences in {1:.2f} seconds".format(len(seqs.keys()),duration),file=sys.stderr)
        return
    


if __name__ == '__main__':
    main()
