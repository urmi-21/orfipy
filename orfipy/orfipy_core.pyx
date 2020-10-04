#!python
#cython: language_level=3
"""
Created on Thu Aug 13 20:01:56 2020

@author: usingh
"""
import re
from collections import deque
from operator import itemgetter
#from ctypes import *
import time
import ahocorasick

   
    
cdef struct ORF:
    int start_index
    int stop_index
    int framenum
    int orf_type
    int length
    #char* start_codon
    #char* stop_codon
    #char* seq
        

cpdef start_search(seq,seq_rc,seqname,minlen,maxlen,strand,starts,stops,table,include_stop,partial3,partial5,find_between_stops,out_opts):
    
    cdef int seq_len=len(seq)
    cdef bint search_rev=False
    cdef bint search_fwd=False
    cdef list orfs_as_struct=[]
    cdef list results=[]
    
    #compile results to return
    bedresults=''
    bed12results=''
    seqresults=['','','']
    
    starts=set(starts)
    stops=set(stops)
    
    if find_between_stops:
        partial3=True
        partial5=True
        findstarts=False
    else:
        findstarts=True
        
    if strand=='b':
        search_fwd=True
        search_rev=True
    elif strand=='f':
        search_fwd=True
        search_rev=False
    elif strand=='r':
        search_fwd=False
        search_rev=True
        
    #init out opts
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
    
    
    #print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXx')
    
    #search fwd strand
    if search_fwd:    
        #search stop positions
        A = ahocorasick.Automaton()
        for idx, key in enumerate(stops):
            A.add_word(key, (idx, key))
        A.make_automaton()
        stop_positions=[end_index-3+1 for end_index, (insert_order, original_value) in A.iter(seq)]
        #pad stop positions to find ORFs without stop in sequence
        #e.g. if seq is AAAATGTTT stops in frame 1:[3] make this--> [-3,3,len(seq)]
        stop_positions=[-3,-2,-1]+stop_positions+[seq_len,seq_len+1,seq_len+2]
        start_positions=None
        if findstarts:
            A = ahocorasick.Automaton()
            for idx, key in enumerate(starts):
                A.add_word(key, (idx, key))
            A.make_automaton()
            start_positions=[end_index-3+1 for end_index, (insert_order, original_value) in A.iter(seq)]
            
        #call search
        orfs_as_struct.extend(get_orfsc(start_positions,
             stop_positions,
             seq_len,
             minlen,
             maxlen,
             partial3, 
             partial5,
             find_between_stops))
        #tlist=[i for i in range(1000000)]
        #testloop(tlist)
        
    if search_rev:
        #search stop positions
        A = ahocorasick.Automaton()
        for idx, key in enumerate(stops):
            A.add_word(key, (idx, key))
        A.make_automaton()
        stop_positions=[end_index-3+1 for end_index, (insert_order, original_value) in A.iter(seq_rc)]
        #pad stop positions to find ORFs without stop in sequence
        #e.g. if seq is AAAATGTTT stops in frame 1:[3] make this--> [-3,3,len(seq)]
        stop_positions=[-3,-2,-1]+stop_positions+[seq_len,seq_len+1,seq_len+2]
        start_positions=None
        if findstarts:
            A = ahocorasick.Automaton()
            for idx, key in enumerate(starts):
                A.add_word(key, (idx, key))
            A.make_automaton()
            start_positions=[end_index-3+1 for end_index, (insert_order, original_value) in A.iter(seq_rc)]
            
            
        #call search
        orfs_as_struct.extend(get_orfsc(start_positions,
             stop_positions,
             seq_len,
             minlen,
             maxlen,
             partial3, 
             partial5,
             find_between_stops,
             True))
    
        
      
    
    
    #print('totallen',len(orfs_as_struct))
       
    #if no output specified only return bed
    if not (bed or bed12 or dna or rna or pep):
        #print('stdout')
        bedresults=orfs_to_bed(orfs_as_struct,seq,seq_rc,seqname,seq_len,include_stop)
        return [bedresults,[],[],[],[]]
    
    
    
    if bed:
        
        bedresults=orfs_to_bed(orfs_as_struct,seq,seq_rc,seqname,seq_len,include_stop)
    if bed12:
        
        bed12results=orfs_to_bed12(orfs_as_struct,seq,seq_rc,seqname,seq_len,include_stop)
    
    #get seq results
    if dna or pep or rna:
        seqresults=orfs_to_seq(orfs_as_struct,seq,seq_rc,seqname,seq_len,include_stop,[dna,rna,pep],table)
    
    #dont change the order of results
    #results.append(bed12results)
    #results.append(bedresults)
    #results.append(dnaresults)
    #results.append(rnaresults)
    #results.append(pepresults)
    
    #dont change the order of retrn list
        
    return [bed12results,bedresults,seqresults[0],seqresults[1],seqresults[2]]
    
    #return [[],[],[],[],[]]
    

'''
cdef testloop(list lint1):
    
    cdef list lint=lint1
    s=time.time()
    for i in range(len(lint)):
        t=lint[i]
        pass
    print('t2',time.time()-s)
    
    
    #this is faster
    s=time.time()
    for i in lint:
        t=i
        pass
    print('t1',time.time()-s)
'''  
    
    

cdef get_orfsc(list start_positions,
             list stop_positions,
             int seq_len,
             int minlen,
             int maxlen, 
             bint partial3=False, 
             bint partial5=False,
             bint find_between_stops=False,
             bint rev_com=False):    
    #set minlen
    if minlen < 3:
        minlen=3   
    
    #create struct
    cdef ORF thisORF
    
    cdef list result=[]
    cdef list start_stop_pairs
    cdef list stops_by_frame=None
    cdef list starts_by_frame=None
    cdef int orf_type #0:complete 1: 5pp 2: 3pp
    #cdef bytes this_start_codon
    cdef int this_frame
    cdef int upstream_stop_index
    cdef int current_stop_index
    cdef int current_start_index
    
   
       
    #print('ARGSSS',partial3,partial5,minlen,maxlen)
    
    #start search
    #split by frame
    stops_by_frame=[[],[],[]]
    for stop in stop_positions:
         stops_by_frame[stop%3].append(stop)
    #print('stops_by_frame',stops_by_frame)
    #if search between start and stop
    if not find_between_stops:
        starts_by_frame=[[],[],[]]
        for start in start_positions:
            starts_by_frame[start%3].append(start)
        #print('SBF',starts_by_frame)
    
    #s=time.time()
    #print('start search')
    start_stop_pairs=find_orfs_between_stopsc(stops_by_frame,starts_by_frame,find_between_stops)
    #print('search done:',time.time()-s)
        
    #format results
    #print('len start loop',len(start_stop_pairs))
    
    for pair in start_stop_pairs:
        #orf_type='complete'
        #orf_type=0
        orf_type=pair[2]
        upstream_stop_index=pair[0]
        #if upstream_stop_index < 0:
        #    orf_type=1
        current_start_index=upstream_stop_index+3
        
        
        current_stop_index=pair[1]
        this_frame=current_stop_index%3
        #fix stop index if its out of sequence
        if current_stop_index >= seq_len:
            #orf_type='3-prime-partial'
            #orf_type=2
            #get total len of ORF till end of seq
            total_len_current=seq_len-current_start_index
            #len for multiple of 3
            required_len=total_len_current-total_len_current%3
            current_stop_index=current_start_index+required_len
        
        
        #check length            
        current_length=current_stop_index-current_start_index
        
        if current_length < minlen or current_length > maxlen:
            continue
               
        #add orfs to results
        if (not partial3) and (orf_type == 2):
            continue
        if (not partial5) and (orf_type == 1):
            continue
        
        framenum=this_frame+1
        if rev_com:
            framenum=-1*framenum
        #add ORF details
        #print('otypes',pair[2],orf_type)
        thisORF.start_index=current_start_index
        thisORF.stop_index=current_stop_index
        thisORF.framenum=framenum
        thisORF.orf_type=orf_type
        thisORF.length=current_length
        #thisORF.seq=current_orf_seq
        #print(thisORF)
        #print(thisORF.start_index)
        result.append(thisORF)
    
    #print('res,len,3p,5p',len(result))
    #print(result)
    return (result)
    

cpdef find_orfs_between_stopsc(list stops_by_frame, list starts_by_frame,bint find_between_stops):
    #result will contain pairs of ORFs [start,stop]
    cdef list result=[]
    cdef int upstream_stop
    cdef int this_start
    cdef bint start_found=False
    
    cdef list stop_size=[0,0,0]
    stop_size[0]=len(stops_by_frame[0])-1
    stop_size[1]=len(stops_by_frame[1])-1
    stop_size[2]=len(stops_by_frame[2])-1
    cdef int otype=0
    
    if not find_between_stops:
        #convert to dequeue for faster pop operation
        starts_by_frame[0]=deque(starts_by_frame[0])
        starts_by_frame[1]=deque(starts_by_frame[1])
        starts_by_frame[2]=deque(starts_by_frame[2])
        
    
    
    #for region between stop codons find a start downstream of upstream stop
        
    #for each frame
    for frame in [0,1,2]:
        #print('in frame:',frame)
        stop_positions=stops_by_frame[frame]
        #print('total stops',len(stop_positions))
        for i, current_stop in enumerate(stop_positions):
            #first po is -1,-2 or -3 ignore
            if i == 0:
                continue
            elif i==stop_size[frame]:
                #this is last stop (out of seq) --> all orfs using this are lacking a stop codon
                #print('reached last',otype)
                otype=2
            upstream_stop=stop_positions[i-1]
            #print('usinc XXXXXXXXXX',starts_by_frame)
            if find_between_stops:
                #if searching only by stop codons orfs can be either complete or 3'partial(lacking stop codon)
                #if orf is not 3'partial make it complete
                #print('adding',upstream_stop,current_stop)
                if not otype==2:
                    otype=0
                result.append((upstream_stop,current_stop,otype))
            else:
                #print('starts_by_frame',starts_by_frame)
                #find  start downstream of upstream_stop_position in this_frame
                #last_start_index=last_start_position_index[frame]
                #for i in range(last_start_index+1,len(starts_by_frame[frame])):
                #for i in starts_by_frame[frame]:
                #if this stop is upstream of next available start
                start_found=False
                while starts_by_frame[frame]:
                    if current_stop <= starts_by_frame[frame][0]:
                        break
                    this_start = starts_by_frame[frame].popleft()
                    
                    if  this_start >= upstream_stop:
                        #found a start for current_stop_position
                        #print('fount',upstream_stop,this_start,current_stop)
                        #update upstream stop as start
                        upstream_stop=this_start-3
                        #last_start=this_start
                        #this is complete orf
                        start_found=True
                        #print('Found strt',otype,start_found)
                        break
                
                
                #other wise is stop is in seq but no start otype is 5'partial else complete
                #0 complete; 1: 5 partial; 2: 3 partial; 3: no start no stop
                if (not start_found) and i==stop_size[frame]: #no start no stop
                    #otype=3
                    continue
                elif start_found and i==stop_size[frame]:
                # start but no stop
                    otype=2
                elif not start_found:
                    otype=1
                else:
                    otype=0
                
                result.append((upstream_stop,current_stop,otype))
                #print('adding',upstream_stop,current_stop,otype,start_found)
                #if not otype ==3:
                #    result.append((upstream_stop,current_stop,otype))
            
    #sort results by position
    result=sorted(result, key=itemgetter(0))
    #print('returnlen',len(result))
    return result


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
    


cdef orfs_to_seq(list orfs_struct_list, str seq, str seq_rc, str seq_name, int seqlen, bint include_stop,out_types,dict table):
    cdef ORF orf
    cdef int ostart
    cdef int oend
    cdef list result=[[],[],[]]
    cdef int ind=0
    #out_types is array of bool [0] --> out seq as nuc; 1-->RNA; 2-->peptide
    
    for orf in orfs_struct_list:
        
        #pair is a list [current_start_index,current_stop_index,this_frame,this_start_codon,this_stop_codon,orftype,seqlen]
        #struct is thisORF.start_index,thisORF.stop_index,thisORF.framenum, thisORF.orf_type, thisORF.length
        ind+=1
        
        ostart=orf.start_index
        oend=orf.stop_index
        frame=orf.framenum
        strand='+'
        current_seq=seq
        if frame < 0:
            strand='-'
            current_seq=seq_rc
        
        otype=orf.orf_type
        orf_type=''
        if otype==0:
            orf_type='complete'
        elif otype==1:
            orf_type='5-prime-partial'
        else:
            orf_type='3-prime-partial'
        olen=orf.length
        
        startcodon="NA"
        stopcodon="NA"
        if otype == 0:
            startcodon=current_seq[ostart:ostart+3]
            stopcodon=current_seq[oend:oend+3]
        elif otype==1:
            #5pp->no stop
            stopcodon=current_seq[oend:oend+3]
        else :
            #nostop
            startcodon=current_seq[ostart:ostart+3]
        #if include stop add 3 to stop pos
        if (otype==0 or otype==1) and include_stop:
            oend+=3
        #get seq
        thisseq_dna=current_seq[ostart:oend]
        
        #if rev complement reverse coordinates
        if frame < 0:
            temp=ostart
            ostart=seqlen-oend
            oend=seqlen-temp
                       
        thisorfid=seq_name+"_ORF."+str(ind)+' ['+str(ostart)+'-'+str(oend)+']('+strand+') type:'+orf_type+' length:'+str(olen)+' frame:'+str(frame)+' start:'+startcodon+' stop:'+stopcodon
        
        
        if out_types[0]:
            result[0].append('>'+thisorfid+'\n'+format_fasta(thisseq_dna))
        if out_types[2]:
            #convert to prot
            #for all codons
            thisseq_pep=''
            for s in [thisseq_dna[i: i + 3] for i in range(0, len(thisseq_dna), 3)]:
                try:
                    thisseq_pep+=table[s]
                except KeyError as error:
                   #print("Error unknown codon:"+s+"...translating to x",file=sys.stderr)
                     thisseq_pep+='X' 
            result[2].append('>'+thisorfid+'\n'+format_fasta(thisseq_pep))
        if out_types[1]:
            #conver to RNA
            thisseq_rna=thisseq_dna.replace('T','U')
            result[1].append('>'+thisorfid+'\n'+format_fasta(thisseq_rna))
    
    #print('return',result)
    return ['\n'.join(res) for res in result] #return a list ofstrings
    #return '\n'.join(result)

#compile results in bed
cdef orfs_to_bed(list orfs_struct_list,str seq,str seq_rc,str seq_name, int seqlen, bint include_stop):
    cdef ORF orf
    cdef int ostart
    cdef int oend
    #print(orfs_list)
    result=[]
    ind=0
    
    for orf in orfs_struct_list:
        ind+=1
        #pair is a list [current_start_index,current_stop_index,this_frame,this_start_codon,this_stop_codon,orftype,seqlen]
        #struct is thisORF.start_index,thisORF.stop_index,thisORF.framenum, thisORF.orf_type, thisORF.length
        #print(orf)
        ostart=orf.start_index
        oend=orf.stop_index
        frame=orf.framenum
        strand='+'
        current_seq=seq
        if frame < 0:
            strand='-'
            current_seq=seq_rc
        
        otype=orf.orf_type
        orf_type=''
        if otype==0:
            orf_type='complete'
        elif otype==1:
            orf_type='5-prime-partial'
        else:
            orf_type='3-prime-partial'
        olen=orf.length
        
        startcodon="NA"
        stopcodon="NA"
        if otype == 0:
            startcodon=current_seq[ostart:ostart+3]
            stopcodon=current_seq[oend:oend+3]
        elif otype==1:
            #5pp->no stop
            stopcodon=current_seq[oend:oend+3]
        else :
            #nostop
            startcodon=current_seq[ostart:ostart+3]
        
        #if include stop add 3 to stop pos
        if (otype==0 or otype==1) and include_stop:
            print(oend,oend+3)
            oend+=3
        
        #if rev complement reverse coordinates
        if frame < 0:
            temp=ostart
            ostart=seqlen-oend
            oend=seqlen-temp
        
        
        thisorfid=seq_name+"_ORF."+str(ind)
        oid= 'ID='+thisorfid+';ORF_type='+orf_type+';ORF_len='+str(olen)+';ORF_frame='+str(frame)+';Start:'+startcodon+';Stop:'+stopcodon
        thisorf=seq_name+'\t'+str(ostart)+'\t'+str(oend)+'\t'+oid+'\t'+'0'+'\t'+strand
        #print(thisorf)
        result.append(thisorf)
    #return as string    
    #print('RETBED:','\n'.join(result))
    return '\n'.join(result)
    #return ''

#compile results in bed12
cdef orfs_to_bed12(list orfs_struct_list,str seq,str seq_rc, str seq_name, int seqlen, bint include_stop):
    cdef ORF orf
    cdef int ostart
    cdef int oend
    #print(orfs_list)
    result=[]
    ind=0
    
    for orf in orfs_struct_list:
        ind+=1
        #pair is a list [current_start_index,current_stop_index,this_frame,this_start_codon,this_stop_codon,orftype,seqlen]
        #struct is thisORF.start_index,thisORF.stop_index,thisORF.framenum, thisORF.orf_type, thisORF.length
        #print(orf)
        ostart=orf.start_index
        oend=orf.stop_index
        frame=orf.framenum
        strand='+'
        current_seq=seq
        if frame < 0:
            strand='-'
            current_seq=seq_rc
        
        otype=orf.orf_type
        orf_type=''
        if otype==0:
            orf_type='complete'
        elif otype==1:
            orf_type='5-prime-partial'
        else:
            orf_type='3-prime-partial'
        olen=orf.length
        
        startcodon="NA"
        stopcodon="NA"
        if otype == 0:
            startcodon=current_seq[ostart:ostart+3]
            stopcodon=current_seq[oend:oend+3]
        elif otype==1:
            #5pp->no stop
            stopcodon=current_seq[oend:oend+3]
        else :
            #nostop
            startcodon=current_seq[ostart:ostart+3]
        
        #if include stop add 3 to stop pos
        if (otype==0 or otype==1) and include_stop:
            oend+=3
                      
        #if rev complement reverse coordinates
        if frame < 0:
            temp=ostart
            ostart=seqlen-oend
            oend=seqlen-temp
            
        thisorfid=seq_name+"_ORF."+str(ind)
        oid= 'ID='+thisorfid+';ORF_type='+orf_type+';ORF_len='+str(olen)+';ORF_frame='+str(frame)+';Start:'+startcodon+';Stop:'+stopcodon
        thisorf=seq_name+'\t'+str(0)+'\t'+str(seqlen)+'\t'+oid+'\t'+'0'+'\t'+strand+'\t'+str(ostart)+'\t'+str(oend)+'\t'+'0'+'\t'+'1'+'\t'+str(seqlen)+'\t'+str(0)
        #print(thisorf)
        result.append(thisorf)
    #return as string    
    return '\n'.join(result)
        
