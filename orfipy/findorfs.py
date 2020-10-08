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
#import pyximport; pyximport.install()
import orfipy_core as oc
from pyfaidx import Fasta
import subprocess

   
    
def worker_map(arglist):
    """
    start worker
    """
    #call orf function
    #poolargs contains [thisseq,thisseq_rc,thisname,minlen,maxlen,strand,starts,stops,include_stop, partial3, partial5, outputs,tmpdir]
    res=oc.start_search(*arglist[:-1])
    
    #directly write res to files
    
    #first open file streams to write the results
    #arglist[-2] is the list of output types
    #for all true values in arglist[-2] open a files in tmp dir and write results
    file_streams=()
    
    for x in range(len(arglist[-2])):
        if arglist[-2][x]:
            #open stream in tmpdir
            file_streams+=(open(os.path.join(arglist[-1],arglist[2]+'.orfipytmp_'+str(x)),'w'),)
        else:
            file_streams+=(None,)
    
    #write results to tmp dir
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
    #poolargs contains [thisseq,thisseq_rc,thisname,minlen,maxlen,strand,starts,stops, table, include_stop, partial3, partial5, outputs,tmpdir]
    #call orf function
    #res=oc.start_search(arglist[0],arglist[1],arglist[2],arglist[3],arglist[4],arglist[5],arglist[6],arglist[7])
    #pass all but last argument
    res=oc.start_search(*arglist[:-1])
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
   
    #return results
    return results_inner

    

def start_multiprocs(seqs, 
                     minlen, 
                     maxlen, 
                     procs,
                     chunk_size, 
                     strand, 
                     starts,
                     stops, 
                     table,include_stop,
                     partial3,
                     partial5,
                     bw_stops,
                     file_streams,
                     tmpdir):
    
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
        poolargs.append([thisseq,thisseq_rc,thisname,minlen,maxlen,strand,starts,stops, table, include_stop, partial3, partial5, bw_stops, outputs,tmpdir])
        
        #if total_read_bytes is more than memory limit
        if total_read_bytes+1000000 >= _MEMLIMIT:
            #mem=float(process.memory_info().rss)/1000000.0
            #print('Processed {0:.2f} MB:'.format(cummulative_read_bytes/1000000),'this c',total_read_bytes,'Memory usage: {0:.2f}MB'.format(mem), end="\r", flush=True,file=sys.stderr )
            #print('Processing {0:d} bytes:'.format(cummulative_read_bytes),'Current Memory usage: {0:.2f}MB'.format(mem), end="\r", flush=True,file=sys.stderr )
            print('Processing {0:d} bytes'.format(cummulative_read_bytes), end="\r", flush=True,file=sys.stderr)
            #process seqs that were added to poolargs
            
            
            
            # find ORFs in currently read data
            
            #if num seq in chunk size are less --> larger seqs; call star_map
            if len(poolargs) < procs-2:
                #print('starting map')
                #results are written to temp files by each worker
                start_map(poolargs,procs)
            else:            
                #call imap unorderd for multiple smaller seqs
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
        print('Processing {0:d} bytes'.format(cummulative_read_bytes), end="\r", flush=True,file=sys.stderr)
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
    

def worker_single(seqs,minlen,maxlen,strand,starts,stops,table,include_stop,partial3,partial5,bw_stops,file_streams,tmp):
    """
    Perform sequential processing

    Parameters
    ----------
    infasta : TYPE
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
    #outputs
    outputs=[]
    for f in file_streams:
        if f:
            outputs.append(True)
        else:
            outputs.append(False)
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
        
        res=oc.start_search(thisseq,thisseq_rc,thisname,minlen,maxlen,strand,starts,stops,table, include_stop, partial3, partial5, bw_stops,outputs)
        
        write_results_single(res, file_streams)


def write_results_single(results,file_streams):
    #results is a list on n lists. each n lists contain a string; file_streams contain n streams to write n lists
    
    all_none=True
    for i in range(len(file_streams)):
        if file_streams[i]:
            if results[i]:
                file_streams[i].write(results[i]+'\n')
            all_none=False
    #no output file specified, print bed results
    if all_none and results[0]:
        print(results[0])

def write_results_multiple(results,file_streams):
    #results is a generator; each item in result is a list of n lists.
    #file_streams contain n file_streams for each list in a result item
    #print(results)
    for reslist in results:
        write_results_single(reslist, file_streams)
    return
    
    
def init_result_files(fileslist,tmp=""):
    #create outdir
    create_outdir(fileslist,tmp)
    #create empty files to append later
    fstreams=()
    for f in fileslist:
        if f:
            fstreams=fstreams+(open(os.path.join(tmp, f),'w'),)
        else:
            fstreams=fstreams+(None,)
    return fstreams


def create_outdir(outlist,outdir):
    for f in outlist:
        if f:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            return

def close_result_files(fstreams):
    for f in fstreams:
        if f:
            f.close()
    
def concat_resultfiles(fstreams,outdir):
    """
    Merge any temporary files, if created. uses subprocess and shell
    """
    #os.chdir(outdir)    
    for f in fstreams:
        if f:
            thisfilename=f.name
            x=fstreams.index(f)
            
            cmd='cat '+outdir+'/*.orfipytmp_'+str(x)+' >> '+thisfilename
            proc = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err = proc.communicate()
            
            cmd='rm '+outdir+'/*.orfipytmp_'+str(x)
            proc = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err = proc.communicate()


def filter_bed_longest(bedfile):
    res={}
    with open(bedfile) as f:
        for ind, l in enumerate(f):
            temp=l.split('\t')
            thisname=temp[0]
            thislen=int(temp[3].split(';')[2].split('=')[-1])
            if thisname in res:
                prevlen=res[thisname][0]
                if thislen > prevlen:
                    res[thisname]=[thislen,l]
            else:
                res[thisname]=[thislen,l]
    
    #write to file
    outfile=os.path.splitext(bedfile)[0]+"_longest"+os.path.splitext(bedfile)[1]
    f=open(outfile,'w')
    for k in res.keys():
        f.write(res[k][1])
    f.close()

def group_bed_frame(bedfile):
    basename=os.path.splitext(bedfile)[0]
    ext=os.path.splitext(bedfile)[1]
    fp1=open(basename+"_+1"+ext,'w')
    fp2=open(basename+"_+2"+ext,'w')
    fp3=open(basename+"_+3"+ext,'w')
    fp4=open(basename+"_-1"+ext,'w')
    fp5=open(basename+"_-2"+ext,'w')
    fp6=open(basename+"_-3"+ext,'w')
    streams={"1":fp1,"2":fp2,"3":fp3,"-1":fp4,"-2":fp5,"-3":fp6}
    with open(bedfile) as f:
        for ind, l in enumerate(f):
            thisframe=l.split('\t')[3].split(';')[3].split('=')[-1]
            streams[thisframe].write(l)
    for k in streams.keys():
        streams[k].close()
    return streams
        
    
    
def group_by_frame_length(bed,bed12,longest,byframe):
    if longest:
        filter_bed_longest(bed)
    if byframe:
        files_frame=group_bed_frame(bed)
        if longest:
            for k in files_frame.keys():
                filter_bed_longest(files_frame[k].name)
            

    
    
    
##########main################
#TODO: handle longest and byframe opts
def main(infasta,
         minlen,
         maxlen,
         procs,
         single_mode,
         chunk_size,
         strand,
         starts,
         stops,
         table,
         include_stop,
         partial3,
         partial5,
         bw_stops,
         longest,
         byframe,
         bed12,
         bed,
         dna,
         rna,
         pep,
         outdir):
    
    """
    

    Parameters
    ----------
    infasta : string
        input path to input fasta file
    minlen : TYPE
        DESCRIPTION.
    maxlen : TYPE
        DESCRIPTION.
    procs : TYPE
        DESCRIPTION.
    single_mode : TYPE
        DESCRIPTION.
    chunk_size : TYPE
        DESCRIPTION.
    strand : TYPE
        DESCRIPTION.
    starts : TYPE
        DESCRIPTION.
    stops : TYPE
        DESCRIPTION.
    nested : TYPE
        DESCRIPTION.
    partial3 : TYPE
        DESCRIPTION.
    partial5 : TYPE
        DESCRIPTION.
    byframe : TYPE
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
    outdir : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    ##start time
    start = time.time()
    
         
    file_streams=init_result_files((bed12, bed, dna, rna, pep),tmp=outdir)    
    
    if not procs:
        procs=int(multiprocessing.cpu_count()*.7)+1
    
    #estimate chunk_size
    if not chunk_size:
        #total mem in MB
        total_mem_MB=psutil.virtual_memory()[0]/1000000
        chunk_size=int(total_mem_MB/(procs*4))
    else:
        chunk_size=int(chunk_size)
    #limit chunk size to 1000 if sequences are extracted; this works best
    if (dna or rna or pep) and (chunk_size > 1000):
        chunk_size=1000
    
    #check py < 3.8; if yes max chunk size can be 2000 otherwise error is reported
    if sys.version_info[1] < 8 and chunk_size > 2000:
        chunk_size = 1900
                
    print("Setting chunk size {} MB. Procs {}".format(chunk_size,procs),file=sys.stderr)
    
    
    #read fasta file
    seqs = Fasta(infasta)
    
    
    
    if single_mode:
        
        worker_single(seqs,
                      minlen,
                      maxlen,
                      strand,
                      starts,
                      stops,
                      table,
                      include_stop,
                      partial3,
                      partial5,
                      bw_stops,
                      file_streams,
                      outdir)
        duration = time.time() - start
    else:
        start_multiprocs(seqs, 
                         minlen,
                         maxlen,
                         procs,
                         chunk_size,
                         strand,
                         starts,
                         stops,
                         table, 
                         include_stop,
                         partial3,
                         partial5,
                         bw_stops,
                         file_streams,
                         outdir)
        duration = time.time() - start
             
               
    
    close_result_files(file_streams)
    #print("Concat...",file=sys.stderr)
    concat_resultfiles(file_streams,outdir)
    
    #after writing all files, write additional file for longest and bystrand
    if longest or byframe:
        bedfile=os.path.join(outdir,bed)
        bed12file=''
        if bed12:
           bed12file=os.path.join(outdir,bed12)
        group_by_frame_length(bedfile,bed12file,longest,byframe)
        
    print("Processed {0:d} sequences in {1:.2f} seconds".format(len(seqs.keys()),duration),file=sys.stderr)


if __name__ == '__main__':
    main()
