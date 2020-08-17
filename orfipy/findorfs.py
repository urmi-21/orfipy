#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 20:01:56 2020

@author: usingh
"""
import time
import pyximport; pyximport.install()
import orfipy.orfipy_core as oc
import sys
from pyfaidx import Fasta
import multiprocessing
from contextlib import closing
#import gc



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



def list_chunks(inlist, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(inlist), n):
        yield inlist[i:i + n]
        

##########main################
def main(infasta,minlen,procs,strand,starts,stops,bed12,bed,dna,rna,pep):
    print("orfipy")
    ##start time
    start = time.time()
    
    single_mode=False
    if not procs:
        single_mode=True
        #procs=multiprocessing.cpu_count()
    
    #read fasta file
    seqs = Fasta(infasta)
    #totalseqs=len(seqs.keys())
    duration = time.time() - start
    print("read in {0:.2f} seconds".format(duration),file=sys.stderr)
    
    if single_mode:
        print('running single')
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
            
    else:
        results=[]
        #poolargsset=[]
        #split data for mp
        poolargs=[]
        #loads all data at once, this can cause problems with big files
        #FIX this
        _LIMIT=50000 #1 GB
        total_read=0
        for s in list(seqs.keys()):
            thisname=s
            thisseq=str(seqs[s])
            this_read=len(thisseq.encode('utf-8'))
            #ignore if seq is < minlen
            if len(thisseq)<minlen:
                continue
            thisseq_rc=None
            if strand == 'b' or strand =='r':
                thisseq_rc=seqs[s][:].complement.reverse.seq
                this_read+=len(thisseq_rc.encode('utf-8'))
                
            
            #check how much is read
            
            if total_read+this_read > _LIMIT:
                print('Running',total_read,'bytes'+'total',total_read+this_read, end="\r", flush=True )
                #print(time.ctime(), end="\r", flush=True)
                #sys.stdout.write("\r" + time.ctime())
                #process current seqs
                with closing( multiprocessing.Pool(procs) ) as p:
                    results_inner = p.imap_unordered(worker2, poolargs, 100)
                for r in results_inner:
                    results.append(r)
                print('Finish Running')
                #poolargsset.append(poolargs)
                del poolargs
                poolargs=[]
                total_read=this_read
            poolargs.append([thisseq,thisseq_rc,thisname,minlen,strand,starts,stops,bed12,bed,dna,rna,pep])
            
        #afterloop process remaining
        with closing( multiprocessing.Pool(procs) ) as p:
            results_inner = p.imap_unordered(worker2, poolargs, 100)
            for r in results_inner:
                results.append(r)
        #poolargsset.append(poolargs)
        
        #print(results)
        
        print('XXXXXXXXXXX')
        #print(len(poolargsset))
        duration = time.time() - start
        print("split in {0:.2f} seconds".format(duration),file=sys.stderr)
    
    
        #with closing( multiprocessing.Pool(procs) ) as p:
        #    results = p.imap_unordered(worker2, poolargs, 100)
        #i=0
        #for poolargs in poolargsset:
        #    i+=1
        #    print('processing',i)
        #    with closing( multiprocessing.Pool(procs) ) as p:
        #        results = p.imap_unordered(worker2, poolargs, 100)
                
        
    #print('total res',(results)) #should be equal to procs
    #for i in results:
    #    print (i)
    #print(results[0])
    start2=time.time()
    
    
    #parse results
    if not (bed12 or bed or dna or rna or pep):
        #print('stdout')
        for reslist in results:
            print(reslist[0]) #reslist[0] contains bed
    else:
        if bed:
            bf=open(bed,'w')
        if bed12:
            b12f=open(bed12,'w')
        if dna:
            dnaf=open(dna,'w')
        if rna:
            rnaf=open(rna,'w')
        if pep:
            pepf=open(pep,'w')
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
        
        
    
    duration = time.time() - start2
    print("write in {0:.2f} seconds".format(duration),file=sys.stderr)
    
    duration = time.time() - start
    print("Processed {0:d} sequences in {1:.2f} seconds".format(len(seqs.keys()),duration),file=sys.stderr)
    
    
    
    

if __name__ == '__main__':
    main()
