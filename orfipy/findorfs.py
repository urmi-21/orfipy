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


'''
def start_workers(poolargs,procs):
    
    #if large seqs start map
    if len(poolargs) < procs-2:
        #print('starting map')
        start_map(poolargs,procs)
        
    else:            
        #print('starting Imap')
        results=start_imap_unordered(poolargs, procs)
'''       
    
def worker_map(arglist):
    """
    start worker
    """
    #call orf function
    #poolargs contains [thisseq,thisseq_rc,thisname,minlen,strand,starts,stops,outputs,tmpdir]
    res=oc.find_orfs(arglist[0],arglist[1],arglist[2],arglist[3],arglist[4],arglist[5],arglist[6],arglist[7])
    
    #directly write to files
    
    file_streams=()
    for x in range(len(arglist[7])):
        if arglist[7][x]:
            #open stream in tmp dir
            file_streams+=(open(os.path.join(arglist[-1],arglist[2]+'.orfipytmp_'+str(x)),'w'),)
        else:
            file_streams+=(None,)
    
    #write_map_results(res,bed12,bed,dna,rna,pep,arglist[12])
    write_results_single(res,file_streams)
    
    del res
    gc.collect()


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
    #poolargs contains [thisseq,thisseq_rc,thisname,minlen,strand,starts,stops,outputs,tmpdir]
    #call orf function
    res=oc.find_orfs(arglist[0],arglist[1],arglist[2],arglist[3],arglist[4],arglist[5],arglist[6],arglist[7])
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
    #bed12=poolargs[0][7]
    #bed=poolargs[0][8]
    #dna=poolargs[0][9]
    #rna=poolargs[0][10]
    #pep=poolargs[0][11]
    #tmpdir=poolargs[0][12]
    #write_results_multiple(results_inner,bed12,bed,dna,rna,pep,tmpdir)
    #return results
    return results_inner

    

def start_multiprocs(seqs,minlen,procs,chunk_size,strand,starts,stops,file_streams,tmpdir):
    
    #outputs
    outputs=[]
    for f in file_streams:
        if f:
            outputs.append(True)
        else:
            outputs.append(False)
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
        poolargs.append([thisseq,thisseq_rc,thisname,minlen,strand,starts,stops,outputs,tmpdir])
        
        #if total_read_bytes is more than memory limit
        if total_read_bytes+1000000 >= _MEMLIMIT:
            #mem=float(process.memory_info().rss)/1000000.0
            #print('Processed {0:.2f} MB:'.format(cummulative_read_bytes/1000000),'this c',total_read_bytes,'Memory usage: {0:.2f}MB'.format(mem), end="\r", flush=True,file=sys.stderr )
            #print('Processing {0:d} bytes:'.format(cummulative_read_bytes),'Current Memory usage: {0:.2f}MB'.format(mem), end="\r", flush=True,file=sys.stderr )
            print('Processed {0:d} bytes'.format(cummulative_read_bytes), end="\r", flush=True,file=sys.stderr)
            #process seqs that were added to poolargs
            
            
            #start_workers(poolargs,procs)
            # find ORFs in currently read data
            if len(poolargs) < procs-2:
                #print('starting map')
                #results are written to temp files by each worker
                start_map(poolargs,procs)
            else:            
                #print('starting Imap')
                results=start_imap_unordered(poolargs, procs)
                #collect and write these results
                write_results_multiple(results,file_streams)
                
            
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
        if len(poolargs) < procs-2:
            #print('starting map')
            #results are written to temp files by each worker
            start_map(poolargs,procs)
        else:            
            #print('starting Imap')
            results=start_imap_unordered(poolargs, procs)
            #collect and write these results
            write_results_multiple(results,file_streams)
        #perform GC
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


def write_results_single(results,file_streams):
    #results is a list on n lists. each n lists contain a string; file_streams contain n streams to write n lists
    
    all_none=True
    for i in range(len(file_streams)):
        if file_streams[i]:            
            file_streams[i].write(results[i]+'\n')
            all_none=False
    #no output file specified, print bed results
    if all_none:
        print(results[0])

def write_results_multiple(results,file_streams):
    #results is a generator; each item in result is a list of n lists.
    #file_streams contain n file_streams for each list in a result item
    for reslist in results:
        write_results_single(reslist, file_streams)
    return
    
    
def init_result_files(fileslist,tmp=""):
    #create empty files to append later
    fstreams=()
    for f in fileslist:
        if f:
            fstreams=fstreams+(open(os.path.join(tmp, f),'w'),)
        else:
            fstreams=fstreams+(None,)
    return fstreams

def close_result_files(fstreams):
    for f in fstreams:
        if f:
            f.close()
    print('closed all')
    
def concat_resultfiles(fstreams,outdir):
    """
    Merge any temporary files, if created
    """
    #os.chdir(outdir)    
    for f in fstreams:
        if f:
            thisfilename=f.name
            x=fstreams.index(f)
            
            cmd='cat '+outdir+'/*.orfipytmp_'+str(x)+' >> '+thisfilename
            print('now',cmd)
            proc = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err = proc.communicate()
            cmd='rm '+outdir+'/*.orfipytmp_'+str(x)
            print('now',cmd)
            proc = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err = proc.communicate()
    
            
    
##########main################
def main(infasta,minlen,procs,single_mode,chunk_size,strand,starts,stops,bed12,bed,dna,rna,pep,outdir):
    ##start time
    start = time.time()
    
    if not outdir:
        outdir="orfipy_"+os.path.basename(infasta)+'_out'
        #print("Temp dir is {}".format(tmpdir),file=sys.stderr)
    
    #create the tmpdir; all tmp out will be stored here
    if not os.path.exists(outdir):
        os.makedirs(outdir)
      
    file_streams=init_result_files((bed12, bed, dna, rna, pep),tmp=outdir)    
    
    
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
    
    
    #read fasta file
    seqs = Fasta(infasta)
    #totalseqs=len(seqs.keys())
    
    
    if single_mode:
        worker_single(seqs, minlen, procs, strand, starts, stops, bed12, bed, dna, rna, pep,outdir)
        duration = time.time() - start
    else:
        start_multiprocs(seqs, minlen, procs, chunk_size, strand, starts, stops, file_streams, outdir)
        duration = time.time() - start
             
               
    
    close_result_files(file_streams)
    print("Concat...",file=sys.stderr)
    concat_resultfiles(file_streams,outdir)
        
    print("Processed {0:d} sequences in {1:.2f} seconds".format(len(seqs.keys()),duration),file=sys.stderr)


if __name__ == '__main__':
    main()
