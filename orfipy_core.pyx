import re
import time

def find_orfs(seq,seq_rc,seqname,minlen=30):
    fwd_res=get_orfs(seq,seqname)
    rev_res=get_orfs(seq_rc,seqname,True)
    combined=fwd_res+rev_res

    #get result as fasta and bed
    results=result_to_bed({seqname:combined},seq,seq_rc)
    print(results)
    #print('XXXXXXXXXXX')
    #print(results[1])
    results=combined
    return results

def find_orfs2(seq,seq_rc,seqname,minlen=30):
    """
    min length is excluding stop
    """
    fwd_res=get_orfs_new(seq,seqname,minlen)
    rev_res=get_orfs_new(seq_rc,seqname,minlen,rev_com=True)
    combined_orfs=fwd_res[0]+rev_res[0]
    combined_seq=fwd_res[1]+rev_res[1]

    #get result as fasta and bed
    results=combined_orfs
    #print('\n'.join(combined_seq))
    #print('XXXXXXXXXXX')
    #print(results[1])
    orfs_to_bed12(combined_orfs,seqname,len(seq))
    orfs_to_seq(combined_orfs,combined_seq,seqname)
    return results


def orfs_to_seq(orfs_list,seq_list,seq_name):
    if not len(orfs_list) == len(seq_list):
        print ("Error")
        return None
    ind=0
    for i in range(len(orfs_list)):
        pair=orfs_list[i]
        ind+=1
        #print(pair[0],pair[1])
        ostart=pair[0]
        oend=pair[1]
        frame=pair[-1]+1
        #determine strand
        strand='+'
        if ostart > oend and oend >= 0:
            #reverse frame
            strand='-'
            #swap
            temp=ostart
            ostart=oend
            oend=temp
        if oend == -9:
            strand='-'
            oend='NA'
        elif oend==-1:
            strand='+'
            oend='NA'
        #compute len
        if oend == 'NA':
            olen='NA'
            otype='incomplete'
        else:
            olen=oend-ostart+1
            otype='complete'
                
        thisorfid=seq_name+"_ORF."+str(ind)+' ['+str(ostart)+'-'+str(oend)+']('+strand+') type:'+otype+' length:'+str(olen)+' frame:'+str(frame)
        print('>'+thisorfid+'\n'+seq_list[i])
    

#compile results in bed12
def orfs_to_bed12(orfs_list,seq_name,seqlen):
    print(orfs_list)
    ind=0
    for pair in orfs_list:
        ind+=1
        #print(pair[0],pair[1])
        ostart=pair[0]
        oend=pair[1]
        frame=pair[-1]+1
        
        
        #determine strand
        strand='+'
        if ostart > oend and oend >= 0:
            #reverse frame
            strand='-'
            #swap
            temp=ostart
            ostart=oend
            oend=temp
        if oend == -9:
            strand='-'
            oend='NA'
        elif oend==-1:
            strand='+'
            oend='NA'
        #compute len
        if oend == 'NA':
            olen='NA'
            otype='incomplete'
        else:
            olen=oend-ostart+1
            otype='complete'
        
        
        thisorfid=seq_name+"_ORF."+str(ind)
        oid= 'ID='+thisorfid+';ORF_type='+otype+';ORF_len='+str(olen)+';ORF_frame='+str(frame)
        thisorf=seq_name+'\t'+str(0)+'\t'+str(seqlen-1)+'\t'+oid+'\t'+strand+'\t'+str(ostart)+'\t'+str(oend)+'\t'+str(0)+'\t'+str(seqlen-1)+'\t'+str(0)
        print(thisorf)
        
            
        
    

def get_orfs_new(seq,seqname,minlen,start_codons=['ATG'],stop_codons=['TAA','TAG','TGA'],report_incomplete=True,rev_com=False):
    cdef int seq_len=len(seq)   
    #get start and stop positions
    cdef list start_positions=[]
    cdef list stop_positions=[]
    
    #re is extremely fast
    for c in start_codons:
        start_positions.extend([m.start() for m in re.finditer(c,seq)])
    for c in stop_codons:
        stop_positions.extend([m.start() for m in re.finditer(c,seq)])
    
    #sort start and stops
    start_positions=sorted(start_positions)
    stop_positions=sorted(stop_positions)
    #divide stops by frames
    cdef list stops_by_frame=[[],[],[]]
    cdef int this_frame
    for i in range(len(stop_positions)):
        this_frame=stop_positions[i]%3
        stops_by_frame[this_frame].append(stop_positions[i])
    #print('SP',start_positions)
    #print('STF',stops_by_frame)
    
    #will contain orfs[start,end,frame]
    cdef list complete_orfs=[]
    cdef list incomplete_orfs=[]
    cdef list complete_orfs_seq=[]
    cdef list incomplete_orfs_seq=[]
    #find all orfs
    cdef int starts_found[3]
    starts_found[:]=[-1,-1,-1]#indices starts found in frame 1 2 3
    cdef int last_stop[3]
    last_stop[:]=[-1,-1,-1]
    cdef int current_start_index
    cdef int current_stop_index
    cdef int last_stop_index
    cdef int current_length
    for i in range(len(start_positions)):
        current_start_index=start_positions[i]
        this_frame=current_start_index%3
        last_stop_index=last_stop[this_frame]
        current_stop_found=False
        #print('lsi',last_stop_index)
        
        #if current_start is contained in prev ORF then ignore: means this start codon lies in another ORF in same frame
        if last_stop_index >=0 and current_start_index < stops_by_frame[this_frame][last_stop_index]:
            #print('inside',current_start_index)
            continue
        
                
        for j in range(last_stop_index+1,len(stops_by_frame[this_frame])):
            #print('running',current_start_index)
            current_stop_index=stops_by_frame[this_frame][j]
            last_stop[this_frame]=j
            if current_stop_index > current_start_index:
                #found ORF!! in this_frame
                current_length=current_stop_index+3-current_start_index
                #print('CL',current_length)
                if current_length >= minlen:
                    complete_orfs.append([current_start_index,current_stop_index+2,this_frame]) #+3 to complete ORF n 0 based corrdinate
                    current_orf_seq = seq[current_start_index:current_stop_index+3]
                    complete_orfs_seq.append(current_orf_seq)
                    current_stop_found=True
                    break
            
            
        #failed to find a stop after searching list of stops
        #if no stops left: ORF has start codon but no downstream in-frame stop codon
        if not current_stop_found:
            #print('No stops for',current_start_index)
            if report_incomplete:
                current_length=seq_len-current_start_index+1
                if current_length >= minlen:
                    incomplete_orfs.append([current_start_index,-1,this_frame])
                    current_orf_seq = seq[current_start_index:1+seq_len-current_length%3]
                    incomplete_orfs_seq.append(current_orf_seq)
                
    
    #print('total orfs',len(complete_orfs))
    #print('total orfs',(complete_orfs))
    #print('ic orfs',incomplete_orfs)
    allorfs=complete_orfs+incomplete_orfs
    allorfs_seq=complete_orfs_seq+incomplete_orfs_seq
    #is seq was reverce complemented, invert the coordinates
    #print('all orfs',allorfs)
    if rev_com:
        allorfs=transform_to_sense(allorfs,seq_len)
        
    #print('all orfs',allorfs)
    return (allorfs,allorfs_seq)
        
        
def transform_to_sense(orfs_list,total_len):
    total_len-=1
    for orf in orfs_list:
        orf[0]=total_len-orf[0]
        if orf[1]>=0:
            orf[1]=total_len-orf[1]
        else:
            orf[1]=-9
    return orfs_list    
    


def get_orfs(seq,seqname,rev_com=False):
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


def result_to_bed(result,seq,seq_rc,minlen=30):
    cdef list bed_result=[]

    for key, value in result.items():
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
                #if negative strand
                if orfstart >= orfend:
                    #swap start and end
                    temp=orfstart
                    orfstart=orfend
                    orfend=temp
                    strand='-'
                #filter by len
                thisorf_len=orfend-orfstart+1
                #skip is short
                if thisorf_len < minlen:
                    continue
                #continue is len passes
                totalorfs+=1
                id='ID='+key+'.'+str(totalorfs)+';ORF len:'+str(thisorf_len)+';ORF frame:'+frame
                score='0'
                thisorf='\t'.join([key,chrstart,chrend,id,score,strand,str(orfstart),str(orfend)])
                bed_result.append(thisorf)

    return '\n'.join(bed_result)


def oldresult_to_bed(result,seq,seq_rc,minlen=30):
    bed_result=[]
    seq_result=[]
    for key, value in result.items():
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
