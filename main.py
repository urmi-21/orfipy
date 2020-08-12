import time
import pyximport; pyximport.install()
import orfs
import sys
from pyfaidx import Fasta
import multiprocessing
import math


def worker(seqlist,result):
    """
    start worker
    """
    #for each fasta find orf
    for k in seqlist:
        thisname=k[0]
        thisseq=k[1]
        thisseq_rc=k[2]
        #add results
        fwd_res=orfs.get_orfs(thisseq)
        rev_res=orfs.get_orfs(thisseq_rc,True)
        combined=fwd_res+rev_res
        result[thisname]=combined
    #print('finished',len(result))


def list_chunks(inlist, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(inlist), n):
        yield inlist[i:i + n]

def result_to_bed(result,fastaob,minlen=30):
    bed_result=[]
    seq_result=[]
    for key, value in result.items():
        seq=fastaob[key][:].seq
        seq_rc=fastaob[key][:].reverse.complement.seq

        chrstart='0'
        chrend=str(len(seq)-1)
        totalorfs=0
        #value is a tuple with three lists
        for i in range(len(value)):
            for orf in value[i]:
                frame=str(i+1)
                if i>=3:
                    frame=str(3-(i+1)) #frames -1,-2,-3
                orfstart=orf[0]
                orfend=orf[1]
                #determine strand
                strand='+'
                #fwd strand seq
                thisseq=seq[orfstart:orfend+1]

                #if negative strand
                if orfstart >= orfend:
                    thisseq=seq_rc[(len(seq)-orfstart):(len(seq)-orfend+1)]
                    #swap start and end
                    temp=orfstart
                    orfstart=orfend
                    orfend=temp
                    strand='-'

                #filter by len
                if len(thisseq)<minlen:
                    continue

                totalorfs+=1
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

def main():
    print("start")
    start = time.time()
    threads=int(sys.argv[2])
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
        print('a')
        thisseqlist=[]
        for k in splitlist[i]:
            thisname=seqs[k].name
            thisseq=seqs[k][:].seq
            thisseq_rc=seqs[k][:].complement.reverse.seq
            thisseqlist.append((thisname,thisseq,thisseq_rc))
        print('b')
        print(i,len(splitlist))
        p = multiprocessing.Process(target=worker, args=(thisseqlist,result))
        jobs.append(p)
        p.start()

    print('wait join')
    for proc in jobs:
        proc.join()
    print('done join')
    duration = time.time() - start
    print('d0',duration)
    #print ('result',result)
    fres=result_to_bed(result,seqs)
    #print (fres[0])
    #print (fres[1])
    #print(result_to_seq(result,seqs))
    duration = time.time() - start
    print('d1',duration)

    #write to file
    bedfile=infasta.split('.')[0]+'.orfs.bed'
    orfsfile=infasta.split('.')[0]+'.orfs.fasta'
    open(bedfile,'w').write(fres[0])
    open(orfsfile,'w').write(fres[1])

    duration = time.time() - start
    print(duration)


if __name__ == '__main__':
    main()
