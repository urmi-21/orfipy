#!python
#cython: language_level=3
"""
Created on Thu Aug 13 20:01:56 2020

@author: usingh
"""
import re
from collections import deque
from operator import itemgetter
from ctypes import *
import time
#from libcpp cimport bool

   
    
cdef struct ORF:
    int start_index
    int stop_index
    int framenum
    int orf_type
    int length
    #char* start_codon
    #char* stop_codon
    #char* seq
        


def get_rev_comp(seq):
    res=seq.replace('A','0').replace('T','A').replace('0','T').replace('G','0').replace('C','G').replace('0','C')[::-1]
    return res

cpdef start_search(seq,
                 seq_rc,
                 seqname,
                 minlen,
                 maxlen,
                 strand,
                 starts,
                 stops,
                 table,
                 include_stop,
                 partial3,
                 partial5, 
                 find_between_stops,
                 out_opts):
    
    
    
    """
    

    Parameters
    ----------
    seq : TYPE
        DESCRIPTION.
    seq_rc : TYPE
        DESCRIPTION.
    seqname : TYPE
        DESCRIPTION.
    minlen : TYPE
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
    TYPE
        Returns a list containing results. each list element is text. if a result type is not specified an empty list is returned at that position
        0 results in bed format
        1 results in bed12 format
        2 dna seq
        3 rna seq
        4 peptide
    """
    cdef int arr[3]
    starts=set(starts)
    stops=set(stops)
    
    
    
    ###Start search
    #re is extremely fast
    '''
    for c in stops:
        stop_positions.extend([m.start() for m in re.finditer(c,seq)])
    #sort start and stops
    stop_positions=sorted(stop_positions)
    #pad stop positions to find ORFs without stop in sequence
    #e.g. if seq is AAAATGTTT stops in frame 1:[3] make this--> [-3,3,len(seq)]
    stop_positions=[-3,-2,-1]+stop_positions+[seq_len,seq_len+1,seq_len+2]
 
    if find_bw_stops:
        #set partial3 5 to true
        partial3=True
        partial5=True    
    elif not find_bw_stops:
        for c in starts:
            start_positions.extend([m.start() for m in re.finditer(c,seq)])
        start_positions=sorted(start_positions)
    '''
    
    
    #print('getting res')
    #out opts
    bed12=out_opts[0]
    bed=out_opts[1]
    dna=out_opts[2]
    rna=out_opts[3]
    pep=out_opts[4]
    #get sequence to write or not
    get_seq=False
    if dna or pep or rna:
        get_seq=True

    if strand=='b':
        #run on fwd strand
        '''
        int *start_positions,
              int *stop_positions,
             int seq_len,
             int minlen,
             int maxlen,
             bint include_stop=False, 
             bint partial3=False, 
             bint partial5=False,
             bint find_bw_stops=False,
             bint return_seqs=False,
             bint rev_com=False
             '''
        
        
        #pyarr = [1, 2, 3, 4]
        #arr = (ctypes.c_int * len(pyarr))(*pyarr)
        
        arr[:]=[1,2,3]
        fwd_res=get_orfs(arr,arr,arr,arr,arr,arr,
                         len(seq),
                         minlen,
                         maxlen,
                         include_stop=include_stop, 
                         partial3=partial3, 
                         partial5=partial5,
                         find_bw_stops=find_between_stops,
                         return_seqs=get_seq)
        #run of rev complemnt strand
        
        #rev_res=get_orfs(seq_rc,seqname,minlen,maxlen,rev_com=True,starts=starts,stops=stops,include_stop=include_stop, partial3=partial3, partial5=partial5,find_bw_stops=find_between_stops,return_seqs=get_seq)
        rev_res=get_orfs(arr,arr,arr,arr,arr,arr,
                         len(seq),
                         minlen,
                         maxlen,
                         include_stop=include_stop, 
                         partial3=partial3, 
                         partial5=partial5,
                         find_bw_stops=find_between_stops,
                         return_seqs=get_seq)
        combined_orfs=fwd_res[0]+rev_res[0]
        combined_seq=fwd_res[1]+rev_res[1]
        
    elif strand == 'f':
        #run on only fwd strand
        #fwd_res=get_orfs(seq,seqname,minlen,maxlen,starts=starts,stops=stops,nested=nested, partial3=partial3, partial5=partial5)
        #fwd_res=get_orfs(seq,seqname,minlen,maxlen,starts=starts,stops=stops,include_stop=include_stop, partial3=partial3, partial5=partial5,find_bw_stops=find_between_stops,return_seqs=get_seq)
        fwd_res=get_orfs(arr,arr,arr,arr,arr,arr,
                         len(seq),
                         minlen,
                         maxlen,
                         include_stop=include_stop, 
                         partial3=partial3, 
                         partial5=partial5,
                         find_bw_stops=find_between_stops,
                         return_seqs=get_seq)
        combined_orfs=fwd_res[0]
        combined_seq=fwd_res[1]
        
    elif strand == 'r':
        #run on only rev comp strand
        #rev_res=get_orfs(seq_rc,seqname,minlen,maxlen,rev_com=True,starts=starts,stops=stops,nested=nested, partial3=partial3, partial5=partial5)
        #rev_res=get_orfs(seq_rc,seqname,minlen,maxlen,rev_com=True,starts=starts,stops=stops,include_stop=include_stop, partial3=partial3, partial5=partial5,find_bw_stops=find_between_stops,return_seqs=get_seq)
        rev_res=get_orfs(arr,arr,arr,arr,arr,arr,
                         len(seq),
                         minlen,
                         maxlen,
                         include_stop=include_stop, 
                         partial3=partial3, 
                         partial5=partial5,
                         find_bw_stops=find_between_stops,
                         return_seqs=get_seq)
        combined_orfs=rev_res[0]
        combined_seq=rev_res[1]
        
    
    
    
    #if no output specified only return bed
    if not (bed or bed12 or dna or rna or pep):
        #print('stdout')
        bedresults=orfs_to_bed(combined_orfs,seqname,len(seq),starts,stops)
        return [bedresults,[],[],[],[]]
    
    #compile results to return
    bedresults=[]
    bed12results=[]
    dnaresults=[]
    rnaresults=[]
    pepresults=[]
    
    if bed:
        bedresults=orfs_to_bed(combined_orfs,seqname,len(seq),starts,stops)
    if bed12:
        bed12results=orfs_to_bed12(combined_orfs,seqname,len(seq),starts,stops)
    if dna:
        dnaresults=orfs_to_seq(combined_orfs,combined_seq,seqname,starts,stops)
    if rna:
        rnaresults=orfs_to_seq(combined_orfs,combined_seq,seqname,starts,stops,out='r')
    if pep:
        pepresults=orfs_to_seq(combined_orfs,combined_seq,seqname,starts,stops,out='p',table=table)
        
    #dont change the order of results
    
    #results.append(bed12results)
    ##results.append(bedresults)
    #results.append(dnaresults)
    #results.append(rnaresults)
    #results.append(pepresults)
    
    
    return [bed12results,bedresults,dnaresults,rnaresults,pepresults]
    
    

#def transcribe_dna(dna):
#    return dna.replace('T','U')

"""
def translate_dna(dna,table):    
    #print(table)
    #use default translation
    if not table:
        #add space after 3 positions; #replace Cys,Gly,Tyr,Ala last
        return " ".join(dna[i: i + 3] for i in range(0, len(dna), 3)).replace('TTT','F').replace('TTC','F').replace('TTA','L').replace('TTG','L').replace('CTT','L').replace('CTC','L').replace('CTA','L').replace('CTG','L').replace('ATT','I').replace('ATC','I').replace('ATA','I').replace('ATG','M').replace('GTT','V').replace('GTC','V').replace('GTA','V').replace('GTG','V').replace('TCT','S').replace('TCC','S').replace('TCA','S').replace('TCG','S').replace('CCT','P').replace('CCC','P').replace('CCA','P').replace('CCG','P').replace('TAT','Y').replace('TAC','Y').replace('TAA','*').replace('TAG','*').replace('CAT','H').replace('CAC','H').replace('CAA','Q').replace('CAG','Q').replace('AAT','N').replace('AAC','N').replace('AAA','K').replace('AAG','K').replace('GAT','D').replace('GAC','D').replace('GAA','E').replace('GAG','E').replace('TGA','*').replace('TGG','W').replace('CGT','R').replace('CGC','R').replace('CGA','R').replace('CGG','R').replace('AGA','R').replace('AGG','R').replace('AGT','S').replace('AGC','S').replace('TGT','0').replace('TGC','0').replace('ACT','1').replace('ACC','1').replace('ACA','1').replace('ACG','1').replace('GCT','2').replace('GCC','2').replace('GCA','2').replace('GCG','2').replace('GGT','3').replace('GGC','3').replace('GGA','3').replace('GGG','3').replace('0','C').replace('1','T').replace('2','A').replace('3','G').replace(' ','')
    
    #if table is provided
    return " ".join(dna[i: i + 3] for i in range(0, len(dna), 3)).replace('TTT',table['TTT']).replace('TTC',table['TTC']).replace('TTA',table['TTA']).replace('TTG',table['TTG']).replace('CTT',table['CTT']).replace('CTC',table['CTC']).replace('CTA',table['CTA']).replace('CTG',table['CTG']).replace('ATT',table['ATT']).replace('ATC',table['ATC']).replace('ATA',table['ATA']).replace('ATG',table['ATG']).replace('GTT',table['GTT']).replace('GTC',table['GTC']).replace('GTA',table['GTA']).replace('GTG',table['GTG']).replace('TCT',table['TCT']).replace('TCC',table['TCC']).replace('TCA',table['TCA']).replace('TCG',table['TCG']).replace('CCT',table['CCT']).replace('CCC',table['CCC']).replace('CCA',table['CCA']).replace('CCG',table['CCG']).replace('TAT',table['TAT']).replace('TAC',table['TAC']).replace('TAA',table['TAA']).replace('TAG',table['TAG']).replace('CAT',table['CAT']).replace('CAC',table['CAC']).replace('CAA',table['CAA']).replace('CAG',table['CAG']).replace('AAT',table['AAT']).replace('AAC',table['AAC']).replace('AAA',table['AAA']).replace('AAG',table['AAG']).replace('GAT',table['GAT']).replace('GAC',table['GAC']).replace('GAA',table['GAA']).replace('GAG',table['GAG']).replace('TGA',table['TGA']).replace('TGG',table['TGG']).replace('CGT',table['CGT']).replace('CGC',table['CGC']).replace('CGA',table['CGA']).replace('CGG',table['CGG']).replace('AGA',table['AGA']).replace('AGG',table['AGG']).replace('AGT',table['AGT']).replace('AGC',table['AGC']).replace('TGT',table['TGT']).replace('TGC',table['TGC']).replace('ACT',table['ACT']).replace('ACC',table['ACC']).replace('ACA',table['ACA']).replace('ACG',table['ACG']).replace('GCT',table['GCT']).replace('GCC',table['GCC']).replace('GCA',table['GCA']).replace('GCG',table['GCG']).replace('GGT',table['GGT']).replace('GGC',table['GGC']).replace('GGA',table['GGA']).replace('GGG',table['GGG']).replace(' ','')
"""    

def format_fasta(seq,int width=62):
    """
    Parameters
    ----------
    seq : TYPE
        DESCRIPTION.
    width : TYPE, num chars on fasta line
        DESCRIPTION. The default is 62. make chars on a line 60

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    return "\n".join(seq[i: i + width] for i in range(0, len(seq), width))
    
    
    
def orfs_to_seq(orfs_list,seq_list,seq_name,starts,stops,out='d',table=None):
    if not len(orfs_list) == len(seq_list):
        print ("Error")
        return ''

    ind=0
    result=[]
    for i in range(len(orfs_list)):
        
        #pair is a list [current_start_index,current_stop_index,this_frame,this_start_codon,this_stop_codon,orftype,seqlen]
        pair=orfs_list[i]
        ind+=1
        #infer strand len etc...
        #ostart,oend,frame,startcodon,stopcodon,strand,otype,olen=format_orf(pair,starts,stops)        
        ostart=pair[0]
        oend=pair[1]
        frame=pair[2]
        strand='+'
        if frame < 0:
            strand='-'
        startcodon=pair[3]
        stopcodon=pair[4]
        otype=pair[5]
        olen=pair[6]
        
        thisorfid=seq_name+"_ORF."+str(ind)+' ['+str(ostart)+'-'+str(oend)+']('+strand+') type:'+otype+' length:'+str(olen)+' frame:'+str(frame)+' start:'+startcodon+' stop:'+stopcodon
        thisseq=seq_list[i]
        if out=='p':
            #convert to prot
            #thisseq=translate_dna(thisseq,table)
            #for all codons
            pepseq=''
            for s in [thisseq[i: i + 3] for i in range(0, len(thisseq), 3)]:
                try:
                    pepseq+=table[s]
                except KeyError as error:
                   #print("Error unknown codon:"+s+"...translating to x",file=sys.stderr)
                     pepseq+='X' 
            thisseq=pepseq

            
        if out == 'r':
            #conver to RNA
            thisseq=thisseq.replace('T','U')
        result.append('>'+thisorfid+'\n'+format_fasta(thisseq))
        #print('>'+thisorfid+'\n'+seq_list[i])
    return '\n'.join(result)

#compile results in bed
def orfs_to_bed(orfs_list,seq_name,seqlen,starts,stops):
    #print(orfs_list)
    result=[]
    ind=0
    for pair in orfs_list:
        ind+=1
        #pair is a list [current_start_index,current_stop_index,this_frame,this_start_codon,this_stop_codon,orftype,seqlen]
        #infer strand len etc...
        #ostart,oend,frame,startcodon,stopcodon,strand,otype,olen=format_orf(pair,starts,stops)
        ostart=pair[0]
        oend=pair[1]
        frame=pair[2]
        strand='+'
        if frame < 0:
            strand='-'
        startcodon=pair[3]
        stopcodon=pair[4]
        otype=pair[5]
        olen=pair[6]
        
        
        
        thisorfid=seq_name+"_ORF."+str(ind)
        oid= 'ID='+thisorfid+';ORF_type='+otype+';ORF_len='+str(olen)+';ORF_frame='+str(frame)+';Start:'+startcodon+';Stop:'+stopcodon
        thisorf=seq_name+'\t'+str(ostart)+'\t'+str(oend)+'\t'+oid+'\t'+'0'+'\t'+strand
        #print(thisorf)
        result.append(thisorf)
    #return as string    
    return '\n'.join(result)

#compile results in bed12
def orfs_to_bed12(orfs_list,seq_name,seqlen,starts,stops):
    #print(orfs_list)
    result=[]
    ind=0
    for pair in orfs_list:
        ind+=1
        #pair is a list [current_start_index,current_stop_index,this_frame,this_start_codon,this_stop_codon,orftype,seqlen]
        #infer strand len etc...
        #ostart,oend,frame,startcodon,stopcodon,strand,otype,olen=format_orf(pair,starts,stops)
        ostart=pair[0]
        oend=pair[1]
        frame=pair[2]
        strand='+'
        if frame < 0:
            strand='-'
        startcodon=pair[3]
        stopcodon=pair[4]
        otype=pair[5]
        olen=pair[6]
        
        
        
        thisorfid=seq_name+"_ORF."+str(ind)
        oid= 'ID='+thisorfid+';ORF_type='+otype+';ORF_len='+str(olen)+';ORF_frame='+str(frame)+';Start:'+startcodon+';Stop:'+stopcodon
        thisorf=seq_name+'\t'+str(0)+'\t'+str(seqlen)+'\t'+oid+'\t'+'0'+'\t'+strand+'\t'+str(ostart)+'\t'+str(oend)+'\t'+'0'+'\t'+'1'+'\t'+str(seqlen)+'\t'+str(0)
        #print(thisorf)
        result.append(thisorf)
    #return as string    
    return '\n'.join(result)
        
#TODO: Write in C
    '''
    int* jagged[2]; 
    // Allocate memory for elements in row 0 
    jagged[0] = malloc(sizeof(int) * 1); 
    // Allocate memory for elements in row 1 
    jagged[1] = malloc(sizeof(int) * 3); 
    '''
cdef get_orfs(int *starts0,int *starts1,int *starts2,
             int *stops0,int *stops1,int *stops2,
             int seq_len,
             int minlen,
             int maxlen,
             bint include_stop=False, 
             bint partial3=False, 
             bint partial5=False,
             bint find_bw_stops=False,
             bint return_seqs=False,
             bint rev_com=False):    
    
    print('XXXXXXXX',starts0[0])
    #set minlen
    if minlen < 3:
        minlen=3
    
    

    
    #create struct
    cdef ORF thisORF
    cdef int orf_type
    #cdef bytes this_start_codon
    
    cdef int this_frame
    
        
        
    if find_bw_stops:
        #report seq between stop codons
        #print('start1')
        #result=find_orfs_between_stops(stops_by_frame)
        #print('end1')
        #set partial3 5 to true
        partial3=True
        partial5=True
        
    else:
        pass
        #result=find_orfs_between_stops(stops_by_frame,starts_by_frame)
                    
    result=[]
    #format results
    #print('start loop')
    for pair in result:
        #orf_type='complete'
        orf_type=0
        upstream_stop_index=pair[0]
        #if upstream_stop_index < 0:
        #    orf_type='5-prime-partial'
        current_start_index=upstream_stop_index+3
        this_start_codon=seq[current_start_index:current_start_index+3]
        
        
        current_stop_index=pair[1]
        this_frame=current_stop_index%3
        #fix stop index if its out of sequence
        if current_stop_index >= seq_len:
            #orf_type='3-prime-partial'
            orf_type=1
            this_stop_codon='NA'
            #get total len of ORF till end of seq
            total_len_current=seq_len-current_start_index
            #len for multiple of 3
            required_len=total_len_current-total_len_current%3
            current_stop_index=current_start_index+required_len
        
        #if start_codon is not in startlist
        if this_start_codon not in starts:
            #orf_type='5-prime-partial'
            orf_type=2
        #get stop codon if present
        if not orf_type == 1:
            this_stop_codon=seq[current_stop_index:current_stop_index+3]
            if len(this_stop_codon)<3:
                this_stop_codon='NA'
        #check length            
        current_length=current_stop_index-current_start_index
        
        if current_length < minlen or current_length > maxlen:
            continue
        #if orf is complete and include_stop
        if include_stop and orf_type == 'complete':
            current_stop_index=current_stop_index+3
        
        #add orfs to results
        if (not partial3) and (orf_type == '3-prime-partial'):
            continue
        if (not partial5) and (orf_type == '5-prime-partial'):
            continue
        #add seq if required
        if return_seqs:
            current_orf_seq = seq[current_start_index:current_stop_index] #current seq
            complete_orfs_seq.append(current_orf_seq)
        framenum=this_frame+1
        if rev_com:
            temp=current_start_index
            current_start_index=seq_len-current_stop_index
            current_stop_index=seq_len-temp
            framenum=-1*framenum
        #add ORF details
        '''
        complete_orfs.append((current_start_index,
                              current_stop_index,
                              framenum,
                              this_start_codon,
                              this_stop_codon,
                              orf_type,
                              current_length))
        '''
        
        thisORF.start_index=current_start_index
        thisORF.stop_index=current_stop_index
        thisORF.framenum=framenum
        thisORF.orf_type=orf_type
        thisORF.length=current_length
        #thisORF.seq=current_orf_seq
        #print(thisORF)
    
        
    return (complete_orfs,complete_orfs_seq)

cdef find_orfs_between_stops(stops_by_frame,starts_by_frame=None):
    if starts_by_frame:
        #conver to dequeue for faster pop operation
        starts_by_frame[0]=deque(starts_by_frame[0])
        starts_by_frame[1]=deque(starts_by_frame[1])
        starts_by_frame[2]=deque(starts_by_frame[2])
    #result will contain pairs of ORFs [start,stop]
    result=[]
    #process for region between stop codons find a start downstream of upstream stop
    last_start_position_index=[-1,-1,-1] #keep track of last start in frame 0,1,2
    last_start=-99999
    
    #process for each frame
    for frame in [0,1,2]:
        #print('in frame:',frame)
        stop_positions=stops_by_frame[frame]
        for i, current_stop in enumerate(stop_positions):
            #first po is -1,-2 or -3 ignore
            if i == 0:
                continue
            #current_stop=stop_positions[i]
            upstream_stop=stop_positions[i-1]
            if not starts_by_frame:
                #print('adding',upstream_stop,current_stop)
                result.append((upstream_stop,current_stop))
            else:
                #find  start downstream of upstream_stop_position in this_frame
                #last_start_index=last_start_position_index[frame]
                #for i in range(last_start_index+1,len(starts_by_frame[frame])):
                #for i in starts_by_frame[frame]:
                #if this stop is upstream of next available start
                while starts_by_frame[frame]:
                    if current_stop <= starts_by_frame[frame][0]:
                        break
                    this_start = starts_by_frame[frame].popleft()
                    
                    if  this_start >= upstream_stop:
                        #found a start for current_stop_position
                        #print('fount',upstream_stop,this_start,current_stop)
                        #update upstream stop as start
                        upstream_stop=this_start-3
                        last_start=this_start
                        break
                #print('adding',upstream_stop,current_stop)
                result.append([upstream_stop,current_stop])
            
    #sort results by position
    result=sorted(result, key=itemgetter(0))
    return result



'''
cdef caller_c(stops_by_frame,starts_by_frame):
    #starts0=NULL
    #starts1=NULL
    #starts2=NULL
    #stops0=NULL
    #stops1=NULL
    #stops2=NULL
    #cdef bint isstart = False
    s=time.time()
    #if starts_by_frame:
    #    isstart=True
    #    starts0 = (c_int * len(starts_by_frame[0]))(*starts_by_frame[0])
    #    starts1 = (c_int * len(starts_by_frame[1]))(*starts_by_frame[1])
    #    starts2 = (c_int * len(starts_by_frame[2]))(*starts_by_frame[2])
    #stops0 = (c_int * len(stops_by_frame[0]))(*stops_by_frame[0])
    #stops1 = (c_int * len(stops_by_frame[1]))(*stops_by_frame[1])
    #stops2 = (c_int * len(stops_by_frame[2]))(*stops_by_frame[2])
    print('Time1xx:',time.time()-s)
    
    #return None
    #return find_orfs_between_stops_c([stops0,stops1,stops2],[starts0,starts1,starts2])
    return find_orfs_between_stops_c(stops_by_frame,starts_by_frame)
'''   
    
#cdef find_orfs_between_stops_c(stops0,stops1,stops2,starts0,starts1,starts2,isstart):
cdef find_orfs_between_stops_c(stops_by_frame,starts_by_frame):
    #s=time.time()
    #result will contain pairs of ORFs [start,stop]
    result=[]
    #process for region between stop codons find a start downstream of upstream stop
    cdef int last_start_position_index[3]
    last_start_position_index[:] = [-1, -1, -1]
    cdef int allframes[3]
    allframes[:]=[0,1,2]
    #last_start_position_index=[-1,-1,-1] #keep track of last start in frame 0,1,2
    cdef int last_start=-99999
    cdef int this_start=-1
    cdef int current_stop=-1
    cdef int upstream_stop=-1
    
    #process for each frame
    for frame in allframes:
        #print('in frame:',frame)
        stop_positions=stops_by_frame[frame]
        for i, current_stop in enumerate(stop_positions):
            #first po is -1,-2 or -3 ignore
            if i == 0:
                continue
            #current_stop=stop_positions[i]
            upstream_stop=stop_positions[i-1]
            if not starts_by_frame:
                #print('adding',upstream_stop,current_stop)
                result.append((upstream_stop,current_stop))
            else:
                #find  start downstream of upstream_stop_position in this_frame
                last_start_index=last_start_position_index[frame]
                for i in range(last_start_index+1,len(starts_by_frame[frame])):
                #for i in starts_by_frame[frame]:
                #if this stop is upstream of next available start
                #while starts_by_frame[frame]:
                    #if current_stop <= starts_by_frame[frame][0]:
                    this_start=starts_by_frame[frame][i]
                    last_start_position_index[frame]+=1
                    if current_stop <= starts_by_frame[frame][i]:
                        break
                    #this_start = starts_by_frame[frame].popleft()
                    
                    if  this_start >= upstream_stop:
                        #found a start for current_stop_position
                        #print('fount',upstream_stop,this_start,current_stop)
                        #update upstream stop as start
                        upstream_stop=this_start-3
                        last_start=this_start
                        break
                #print('adding',upstream_stop,current_stop)
                result.append([upstream_stop,current_stop])
    
    #sort results by position
    result=sorted(result, key=itemgetter(0))
    #print('Time2xx:',time.time()-s)
    return result
    




