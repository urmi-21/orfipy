import sys
import re
import time
from pyfaidx import Fasta
import multiprocessing
import math


def get_orfs(seq):
    start_codons=['ATG']
    stop_codons=['TAA','TAG','TGA']
    #get start and stop positions
    start_positions=[]
    stop_positions=[]

    start = time.time()
    for c in start_codons:
        start_positions.extend([m.start() for m in re.finditer(c,seq)])
    for c in stop_codons:
        stop_positions.extend([m.start() for m in re.finditer(c,seq)])
    duration = time.time() - start
    #print(duration)
    #print(start_positions)
    #print(stop_positions)
    #print('search done',len(start_positions),len(stop_positions))
    #divide into frames
    start_1=[]
    stop_1=[]
    start_2=[]
    stop_2=[]
    start_3=[]
    stop_3=[]
    for p in start_positions:
        x=p%3
        if x==0:
            start_1.append(p)
        elif x==1:
            start_2.append(p)
        elif x==2:
            start_3.append(p)
    for p in stop_positions:
        x=p%3
        if x==0:
            stop_1.append(p)
        elif x==1:
            stop_2.append(p)
        elif x==2:
            stop_3.append(p)
    #print('split done')
    #sort all lists
    start_1=sorted(start_1)
    start_2=sorted(start_2)
    start_3=sorted(start_3)
    stop_1=sorted(stop_1)
    stop_2=sorted(stop_2)
    stop_3=sorted(stop_3)
    #print('sort done')
    #print("by frames")
    #print(start_1,stop_1)
    #print(start_2,stop_2)
    #print(start_3,stop_3)

    #print("start")
    #get orfs for all 3 frames
    orfs_1=get_completeORFs(start_1,stop_1)
    orfs_2=get_completeORFs(start_2,stop_2)
    orfs_3=get_completeORFs(start_3,stop_3)
    return[orfs_1,orfs_2,orfs_3]

def get_completeORFs(starts,stops):
    s=starts
    t=stops
    orfs=[]
    start_index=0
    stop_index=0
    last_stop=-1
    last_start=-1
    last_stop_index=0
    for start_index in range(len(s)):
        #break if no stop codons downstream
        if last_stop_index>=len(t)-1:
            break
        current_start=s[start_index]
        #if this start is contained in an ORF of the same frame then ignore
        if current_start <= last_stop:
            continue
        #search for stop
        for stop_index in range(last_stop_index,len(t)):
            if t[stop_index]>current_start:
                #found stop codon
                last_stop=t[stop_index]
                orfs.append((current_start,last_stop))
                last_stop_index=stop_index
                break
            #no stop found
            last_stop_index=stop_index

    #print(orfs)
    return orfs

def worker(keys,seqs,result):
    """
    start worker
    """
    #for each fasta find orf
    #print(keys)
    for k in keys:
        thisseq=seqs[k][:].seq
        thisname=seqs[k].name
        #add results
        result[thisname]=get_orfs(thisseq)


def list_chunks(inlist, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(inlist), n):
        yield inlist[i:i + n]

def result_to_bed(result,fastaob,minlen=30):
    bed_result=[]
    seq_result=[]
    for key, value in result.items():
        seq=fastaob[key][:]
        totalorfs=0
        #value is a tuple with three lists
        for i in range(len(value)):
            for orf in value[i]:
                frame=str(i+1)
                orfstart=orf[0]
                orfend=orf[1]+2
                #determine strand
                strand='+'
                #if negative strand
                if orfstart >= orfend:
                    strand='-'
                thisseq=seq[orfstart:orfend+1]

                #filter by len
                if len(thisseq)<minlen:
                    continue

                totalorfs+=1
                chrstart='0'
                chrend=str(len(seq)-1)
                id='ID='+key+'.'+str(totalorfs)+';ORF len:'+str(len(thisseq))+';ORF frame:'+frame
                score='0'

                thisorf='\t'.join([key,chrstart,chrend,id,score,strand,str(orfstart),str(orfend)])
                bed_result.append(thisorf)

                ##get seq
                thisname='>'+key+'.'+str(totalorfs)+' '+str(orfstart)+'-'+str(orfend)+'('+strand+')'+' len:'+str(len(thisseq))+' ORF frame:'+frame
                seq_result.append(str(thisname)+'\n'+str(thisseq))
                if not thisseq:
                    print ("EROOORRRRRRRRRRRRRRRRRRR")

    return ('\n'.join(bed_result),'\n'.join(seq_result))

def result_to_seq(result,fastaob,minlen=30):
    seq_result=[]
    for key, value in result.items():
        seq=fastaob[key][:]
        #value is a tuple with three lists
        for i in range(len(value)):
            for orf in value[i]:
                thisname='>'+key+str(orf[0])+'-'+str(orf[1])
                #seq_result.append(thisname)
                thisseq=seq[orf[0]:orf[1]+3]
                if len(thisseq)>=minlen:
                    seq_result.append(str(thisname)+'\n'+str(thisseq))
                if not thisseq:
                    print ("EROOORRRRRRRRRRRRRRRRRRR")

    print(seq_result)
    return '\n'.join(seq_result)


#test="ATACGCATGTTTATGCCATTAGATAGTAAATCAGACTCAGACGACGACTACGACTCAGCATCAGCAGGCAATCAGCACGAGCTACAGCACTCAGACTACGACGACTCAGCAT"
#test='gcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcggcgtaa'.upper()
#with open(sys.argv[1]) as f:
    #d=f.read().splitlines()
#test=''.join(d)
#print(len(test))
#get_orfs(test)
#get_completeORFs([4,12,17,25],[10,13,22,23,33,40])

#read fasta file

#split keys in workers
if __name__ == '__main__':
    threads=4
    infasta=sys.argv[1]
    seqs = Fasta(infasta)
    totalseqs=len(seqs.keys())
    seqperthread=math.ceil((totalseqs/threads))
    splitlist= list(list_chunks(list(seqs.keys()),seqperthread))
    #print('totalseq',totalseqs,'spt',seqperthread,'lst len',len(splitlist))
    #for i in list_chunks(list(seqs.keys()),threads):
    #    splitlist.append(i)
    #print('splitt',splitlist)

    #split keys
    manager = multiprocessing.Manager()
    result = manager.dict()
    jobs = []
    for i in range(len(splitlist)):
        p = multiprocessing.Process(target=worker, args=(splitlist[i],seqs,result))
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()

    #print ('result',result)
    fres=result_to_bed(result,seqs)
    #print(result_to_seq(result,seqs))
