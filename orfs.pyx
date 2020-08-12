import re
import time

def get_orfs(seq,rev_com=False):
    seq_len=len(seq)
    start_codons=['ATG']
    stop_codons=['TAA','TAG','TGA']
    #get start and stop positions
    cdef list start_positions=[]
    cdef list stop_positions=[]
    #start = time.time()
    for c in start_codons:
        start_positions.extend([m.start() for m in re.finditer(c,seq)])
    for c in stop_codons:
        stop_positions.extend([m.start() for m in re.finditer(c,seq)])
    #duration = time.time() - start
    #print(duration)

    #print(start_positions)
    #print(stop_positions)
    #print('search done',len(start_positions),len(stop_positions))
    #divide into frames
    cdef list start_1=[]
    cdef list stop_1=[]
    cdef list start_2=[]
    cdef list stop_2=[]
    cdef list start_3=[]
    cdef list stop_3=[]
    #find translation +1,+2,+3 frame
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

    #if reverse complement, transform coordinates
    if rev_com:
        orfs_1=transform_to_sense(orfs_1,seq_len)
        orfs_2=transform_to_sense(orfs_2,seq_len)
        orfs_3=transform_to_sense(orfs_3,seq_len)

    return[orfs_1,orfs_2,orfs_3]

def transform_to_sense(orfs_list,total_len):
    for orf in orfs_list:
        orf[0]=total_len-orf[0]
        orf[1]=total_len-orf[1]
    return orfs_list


def get_completeORFs(starts,stops):
    cdef list s=starts
    cdef list t=stops
    cdef list orfs=[]
    cdef int start_index=0
    cdef int stop_index=0
    cdef int last_stop=-1
    cdef int last_start=-1
    cdef int last_stop_index=0
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
                #add +2
                orfs.append([current_start,last_stop+2])
                last_stop_index=stop_index
                break
            #no stop found
            last_stop_index=stop_index
      #print(orfs)
    return orfs
