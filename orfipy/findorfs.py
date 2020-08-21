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
import subprocess



def start_workers(poolargs,procs):
    
    #if large seqs start map
    if len(poolargs) < procs-2:
        #print('starting map')
        start_map(poolargs,procs)
        #return res
    else:            
        #print('starting Imap')
        start_imap_unordered(poolargs, procs)
        
    
def worker_map(arglist):
    """
    start worker
    """
    #call orf function
    res=oc.find_orfs(arglist[0],arglist[1],arglist[2],arglist[3],arglist[4],arglist[5],arglist[6],arglist[7],arglist[8],arglist[9],arglist[10],arglist[11])
    
    #print(res)
    
    #poolargs=[thisseq,thisseq_rc,thisname,minlen,strand,starts,stops,bed12,bed,dna,rna,pep,tmpdir]
    #directly write to files
    bed12=arglist[7]
    bed=arglist[8]
    dna=arglist[9]
    rna=arglist[10]
    pep=arglist[11]
    
    if bed12:
        bed12=arglist[2]+'.tmp_bed12'
    if bed:
        bed=arglist[2]+'.tmp_bed'
    if dna:
        dna=arglist[2]+'.tmp_dna'
    if rna:
        rna=arglist[2]+'.tmp_rna'
    if pep:
        pep=arglist[2]+'.tmp_pep'
    
    #write_map_results(res,bed12,bed,dna,rna,pep,arglist[12])
    write_results_single(res,bed12,bed,dna,rna,pep,arglist[12])
    #file=os.path.join(arglist[-1],arglist[2]+"_tmp")
    #f=open(file,'w')
    #f.write(res[2])
    #f.close()
    #del res
    #gc.collect()
    #print('return')
    #return True

#def write_map_results(results,bed12,bed,dna,rna,pep,tmpdir):
#    write_results_single(results,bed12,bed,dna,rna,pep,tmpdir)

def start_map(poolargs,procs):
    """
    Suitable for large sequences with a large number of ORF e.g. genomes
    Low IO overhead

    Parameters
    ----------
    poolargs : TYPE
        DESCRIPTION.
    procs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    pool=multiprocessing.Pool(processes=procs)
    pool.map(worker_map,poolargs)
    pool.close()
    pool.join()
    



def worker_imap(arglist):
    """
    start worker
    """
    #call orf function
    res=oc.find_orfs(arglist[0],arglist[1],arglist[2],arglist[3],arglist[4],arglist[5],arglist[6],arglist[7],arglist[8],arglist[9],arglist[10],arglist[11])
    return res
    
def start_imap_unordered(poolargs,procs):
    """
    Suitable for smaller sequences with a less number of ORF e.g. transcriptomes

    Parameters
    ----------
    poolargs : TYPE
        DESCRIPTION.
    procs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    with closing( multiprocessing.Pool(procs) ) as p:
        results_inner = p.imap_unordered(worker_imap, poolargs, 100)
    
    #return results_inner
    #results inner is generator with results
    #write results to file
    bed12=poolargs[0][7]
    bed=poolargs[0][8]
    dna=poolargs[0][9]
    rna=poolargs[0][10]
    pep=poolargs[0][11]
    tmpdir=poolargs[0][12]
    write_results_multiple(results_inner,bed12,bed,dna,rna,pep,tmpdir)

    

def start_multiprocs(seqs,minlen,procs,chunk_size,strand,starts,stops,bed12,bed,dna,rna,pep,tmpdir):
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
        poolargs.append([thisseq,thisseq_rc,thisname,minlen,strand,starts,stops,bed12,bed,dna,rna,pep,tmpdir])
        
        #if total_read_bytes is more than memory limit
        if total_read_bytes+1000000 >= _MEMLIMIT:
            #mem=float(process.memory_info().rss)/1000000.0
            #print('Processed {0:.2f} MB:'.format(cummulative_read_bytes/1000000),'this c',total_read_bytes,'Memory usage: {0:.2f}MB'.format(mem), end="\r", flush=True,file=sys.stderr )
            #print('Processing {0:d} bytes:'.format(cummulative_read_bytes),'Current Memory usage: {0:.2f}MB'.format(mem), end="\r", flush=True,file=sys.stderr )
            print('Processed {0:d} bytes'.format(cummulative_read_bytes), end="\r", flush=True,file=sys.stderr)
            #process seqs that were added to poolargs
            
            
            start_workers(poolargs,procs)
            
            #perform GC
            #del results_inner
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
        start_workers(poolargs,procs)
        #perform GC
        #del results_inner
        del poolargs
        gc.collect()
        #create empty list
        poolargs=[]
        #reset total read bytes for next batch
        total_read_bytes=0
    
    print()
    

def worker_single(seqs,minlen,procs,strand,starts,stops,bed12,bed,dna,rna,pep,tmp):
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
    #results=[]
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
        write_results_single(res, bed12, bed, dna, rna, pep,tmp)


def write_results_single(results,bed12,bed,dna,rna,pep,tmp=""):
    #print(results)
    #if no out type is specified
    if not (bed12 or bed or dna or rna or pep):
        #print('stdout')
        print(results[0]) #reslist[0] contains bed
    else:
        if bed:
            bf=open(os.path.join(tmp, bed),'a')
        if bed12:
            b12f=open(os.path.join(tmp, bed12),'a')
        if dna:
            dnaf=open(os.path.join(tmp, dna),'a')
        if rna:
            rnaf=open(os.path.join(tmp, rna),'a')
        if pep:
            pepf=open(os.path.join(tmp, pep),'a')
        if bed:
            bf.write(results[0]+'\n')
        if bed12:
            b12f.write(results[1]+'\n')
        if dna:
            dnaf.write(results[2]+'\n')
        if rna:
            rnaf.write(results[3]+'\n')
        if pep:
            pepf.write(results[4]+'\n')
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


def write_results_multiple(results,bed12,bed,dna,rna,pep,tmp=""):
    for reslist in results:
        write_results_single(reslist, bed12, bed, dna, rna, pep,tmp)
    return
    
    
def init_result_files(bed12,bed,dna,rna,pep,tmp=""):
    #create empty files to append later
    if bed:
        f=open(os.path.join(tmp, bed),'w')
        f.close()
    if bed12:
        f=open(os.path.join(tmp, bed12),'w')
        f.close()
    if dna:
        f=open(os.path.join(tmp, dna),'w')
        f.close()
    if rna:
        f=open(os.path.join(tmp, rna),'w')
        f.close()
    if pep:
        f=open(os.path.join(tmp, pep),'w')
        f.close()


def concat_resultfiles(bed12,bed,dna,rna,pep,outdir):
    os.chdir(outdir)
    #allfiles=[file for file in os.listdir('.')]
    #print(allfiles)
    
    #combine bed12
    if bed12:
        proc = subprocess.Popen('cat *tmp_bed12 >> '+bed12, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = proc.communicate()
        proc = subprocess.Popen('rm *tmp_bed12', shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = proc.communicate()
    if bed:
        proc = subprocess.Popen('cat *tmp_bed >> '+bed, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = proc.communicate()
        proc = subprocess.Popen('rm *tmp_bed', shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = proc.communicate()
    if dna:
        proc = subprocess.Popen('cat *tmp_dna >> '+dna, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = proc.communicate()
        proc = subprocess.Popen('rm *tmp_dna', shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = proc.communicate()
    if rna:
        proc = subprocess.Popen('cat *tmp_rna >> '+rna, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = proc.communicate()
        proc = subprocess.Popen('rm *tmp_rna', shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = proc.communicate()
    if pep:
        proc = subprocess.Popen('cat *tmp_pep >> '+pep, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = proc.communicate()
        proc = subprocess.Popen('rm *tmp_pep', shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = proc.communicate()
    
            
    
##########main################
def main(infasta,minlen,procs,single_mode,chunk_size,strand,starts,stops,bed12,bed,dna,rna,pep,outdir):
    
    if not outdir:
        outdir="orfipy_"+os.path.basename(infasta)+'_out'
        #print("Temp dir is {}".format(tmpdir),file=sys.stderr)
    
    #create the tmpdir; all tmp out will be stored here
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    
    init_result_files(bed12, bed, dna, rna, pep,tmp=outdir)    
    
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
    if sys.version_info[1] < 8 and chunk_size > 2000:
        chunk_size = 1900
                
    print("Setting chunk size {} MB".format(chunk_size),file=sys.stderr)
    #sys.exit(1)
    #read fasta file
    seqs = Fasta(infasta)
    #totalseqs=len(seqs.keys())
    
    
    if single_mode:
        worker_single(seqs, minlen, procs, strand, starts, stops, bed12, bed, dna, rna, pep,outdir)
        duration = time.time() - start
        print("Processed {0:d} sequences in {1:.2f} seconds".format(len(seqs.keys()),duration),file=sys.stderr)
        return
    else:
        start_multiprocs(seqs, minlen, procs, chunk_size, strand, starts, stops, bed12, bed, dna, rna, pep,outdir)
        duration = time.time() - start
        print("Processed {0:d} sequences in {1:.2f} seconds".format(len(seqs.keys()),duration),file=sys.stderr)
        
        print("Concat...",file=sys.stderr)
        concat_resultfiles(bed12,bed,dna,rna,pep,outdir)
        return
    


if __name__ == '__main__':
    main()
