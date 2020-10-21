#!python
#cython: language_level=3
"""
Created on Thu Aug 13 20:01:56 2020

@author: usingh
"""
from collections import deque
from operator import itemgetter
import ahocorasick
import itertools

   
cdef struct ORF:
    int start_index
    int stop_index
    int framenum 
    int orf_type
    int length
        

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
    
           
    #if no output specified only return bed
    if not (bed or bed12 or dna or rna or pep):
        bedresults=orfs_to_bed(orfs_as_struct,seq,seq_rc,seqname,seq_len,include_stop)
        return [bedresults,[],[],[],[]]    
    
    if bed:
        
        bedresults=orfs_to_bed(orfs_as_struct,seq,seq_rc,seqname,seq_len,include_stop)
    if bed12:
        
        bed12results=orfs_to_bed12(orfs_as_struct,seq,seq_rc,seqname,seq_len,include_stop)
    
    #get seq results
    if dna or pep or rna:
        seqresults=orfs_to_seq(orfs_as_struct,seq,seq_rc,seqname,seq_len,include_stop,[dna,rna,pep],table)
    
       
    #dont change the order of retrn list
    
    return [bed12results,bedresults,seqresults[0],seqresults[1],seqresults[2]]
    
    
    

    

cdef list get_orfsc(list start_positions,
             list stop_positions,
             int seq_len,
             int minlen,
             int maxlen, 
             bint partial3=False, 
             bint partial5=False,
             bint find_between_stops=False,
             bint rev_com=False):
    """
    """
    
    
    
    
    
    
    #set minlen
    if minlen < 3:
        minlen=3   
    
    #create struct obj
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
    
   
      
    
    #start search
    #split by frame
    stops_by_frame=[[],[],[]]
    for stop in stop_positions:
         stops_by_frame[stop%3].append(stop) 
    
    
    if find_between_stops:
        start_stop_pairs=list(map(find_between_stops_v,stops_by_frame))
    else:
        #if search between start and stop
        starts_by_frame=[[],[],[]]
        for start in start_positions:
            starts_by_frame[start%3].append(start)
            
        start_stop_pairs=list(map(find_between_start_stop_v,starts_by_frame,stops_by_frame))
        
        
    #format results
    for pair in list(itertools.chain.from_iterable(start_stop_pairs)):
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
        thisORF.start_index=current_start_index
        thisORF.stop_index=current_stop_index
        thisORF.framenum=framenum
        thisORF.orf_type=orf_type
        thisORF.length=current_length
    
        result.append(thisORF)   

    return result


    
cpdef list find_between_stops_v(list stop_positions):
    #find regions between stops
    cdef list result=[]
    cdef int upstream_stop
    cdef int current_stop
    cdef int otype=0
    otypes={True:0,False:2}
    cdef int total_stops=len(stop_positions)
    for i in range(1,total_stops):
        upstream_stop=stop_positions[i-1]
        current_stop=stop_positions[i]
        #if i == total_stops-1:
        #    otype=2
        otype=otypes[bool(total_stops-1-i)] #is false at last index
        
        result.append((upstream_stop,current_stop,otype))
    return result
    
cpdef list find_between_start_stop_v(start_positions,  list stop_positions):
    #result will contain pairs of ORFs [start,stop]
    cdef list result=[]
    cdef int upstream_stop
    cdef int this_start
    cdef bint start_found=False
    cdef int otype=0
    cdef int total_stops=len(stop_positions)
    cdef int i=0
    
    #convert to dequeue for faster pop operation
    start_positions=deque(start_positions)
    #print('S',start_positions)
    #print('T',stop_positions)

    
    #first stop pos is -1,-2 or -3 ignore
    for i in range(1,total_stops):
        upstream_stop=stop_positions[i-1]
        current_stop=stop_positions[i]
        #if this stop is upstream of next available start
        start_found=False
        while start_positions:
            if current_stop <= start_positions[0]:
                break
            this_start = start_positions.popleft()
            
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
        if (not start_found) and i==total_stops-1: #no start no stop
            #otype=3
            continue
        elif start_found and i==total_stops-1:
         #start but no stop
            otype=2
        elif not start_found:
            otype=1
        else:
            otype=0
        
        result.append((upstream_stop,current_stop,otype))
    
    return result
    
    


cdef str format_fasta(seq):
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
    cdef int width=62
    return "\n".join(seq[i: i + width] for i in range(0, len(seq), width))
    


cdef list orfs_to_seq(list orfs_struct_list, str seq, str seq_rc, str seq_name, int seqlen, bint include_stop,out_types,dict table):
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

    


#compile results in bed
cdef str orfs_to_bed(list orfs_struct_list,str seq,str seq_rc,str seq_name, int seqlen, bint include_stop):
    cdef ORF orf
    cdef int ostart
    cdef int oend
    
    result=[]
    ind=0
    
    for orf in orfs_struct_list:
        ind+=1
        #pair is a list [current_start_index,current_stop_index,this_frame,this_start_codon,this_stop_codon,orftype,seqlen]
        #struct is thisORF.start_index,thisORF.stop_index,thisORF.framenum, thisORF.orf_type, thisORF.length
        
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
            #print(oend,oend+3)
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
    return '\n'.join(result)
    

#compile results in bed12
cdef str orfs_to_bed12(list orfs_struct_list,str seq,str seq_rc, str seq_name, int seqlen, bint include_stop):
    cdef ORF orf
    cdef int ostart
    cdef int oend

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
        
