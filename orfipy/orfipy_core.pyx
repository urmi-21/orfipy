#!python
#cython: language_level=3
"""
Created on Thu Aug 13 20:01:56 2020

@author: usingh
"""
import re
import time


def get_rev_comp(seq):
    res=seq.replace('A','0').replace('T','A').replace('0','T').replace('G','0').replace('C','G').replace('0','C')[::-1]
    return res

def start_search(seq,seq_rc,seqname,minlen,maxlen,strand,starts,stops,nested, partial3, partial5, out_opts):
    
    print('find_orfs params:',seqname,minlen,maxlen,strand,starts,stops,nested, partial3, partial5, out_opts)
    
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
    #print('getting res')

    if strand=='b':
        #run on fwd strand    
        fwd_res=get_orfs(seq,seqname,minlen,maxlen,starts=starts,stops=stops,nested=nested, partial3=partial3, partial5=partial5)
        #run of rev complemnt strand
        rev_res=get_orfs(seq_rc,seqname,minlen,maxlen,rev_com=True,starts=starts,stops=stops,nested=nested, partial3=partial3, partial5=partial5)
        combined_orfs=fwd_res[0]+rev_res[0]
        combined_seq=fwd_res[1]+rev_res[1]
        
    elif strand == 'f':
        #run on only fwd strand
        fwd_res=get_orfs(seq,seqname,minlen,maxlen,starts=starts,stops=stops,nested=nested, partial3=partial3, partial5=partial5)
        combined_orfs=fwd_res[0]
        combined_seq=fwd_res[1]
        
    elif strand == 'r':
        #run on only rev comp strand
        rev_res=get_orfs(seq_rc,seqname,minlen,maxlen,rev_com=True,starts=starts,stops=stops,nested=nested, partial3=partial3, partial5=partial5)
        combined_orfs=rev_res[0]
        combined_seq=rev_res[1]
        
    
    #out opts
    
    bed=out_opts[0]
    bed12=out_opts[1]
    dna=out_opts[2]
    rna=out_opts[3]
    pep=out_opts[4]
    
    #if no output specified only return bed
    if not (bed or bed12 or dna or rna or pep):
        #print('stdout')
        bedresults=orfs_to_bed12(combined_orfs,seqname,len(seq))
        return [bedresults,[],[],[],[]]
    
    #compile results to return
    results=[]
    
    bedresults=[]
    bed12results=[]
    dnaresults=[]
    rnaresults=[]
    pepresults=[]
    
    if bed:
        bedresults=orfs_to_bed12(combined_orfs,seqname,len(seq))
    if bed12:
        bed12results=orfs_to_bed12(combined_orfs,seqname,len(seq))
    if dna:
        dnaresults=orfs_to_seq(combined_orfs,combined_seq,seqname)
    if rna:
        rnaresults=orfs_to_seq(combined_orfs,combined_seq,seqname,out='r')
    if pep:
        pepresults=orfs_to_seq(combined_orfs,combined_seq,seqname,out='p')
        
    
    results.append(bedresults)
    results.append(bed12results)
    results.append(dnaresults)
    results.append(rnaresults)
    results.append(pepresults)
    
    
    return results
    
    

def transcribe_dna(dna):
    return dna.replace('T','U')

def translate_dna(dna):    
    #add space after 3 positions; #replace Cys,Gly,Tyr,Ala last
    return " ".join(dna[i: i + 3] for i in range(0, len(dna), 3)).replace('TTT','F').replace('TTC','F').replace('TTA','L').replace('TTG','L').replace('CTT','L').replace('CTC','L').replace('CTA','L').replace('CTG','L').replace('ATT','I').replace('ATC','I').replace('ATA','I').replace('ATG','M').replace('GTT','V').replace('GTC','V').replace('GTA','V').replace('GTG','V').replace('TCT','S').replace('TCC','S').replace('TCA','S').replace('TCG','S').replace('CCT','P').replace('CCC','P').replace('CCA','P').replace('CCG','P').replace('TAT','Y').replace('TAC','Y').replace('TAA','.').replace('TAG','.').replace('CAT','H').replace('CAC','H').replace('CAA','Q').replace('CAG','Q').replace('AAT','N').replace('AAC','N').replace('AAA','K').replace('AAG','K').replace('GAT','D').replace('GAC','D').replace('GAA','E').replace('GAG','E').replace('TGA','.').replace('TGG','W').replace('CGT','R').replace('CGC','R').replace('CGA','R').replace('CGG','R').replace('AGA','R').replace('AGG','R').replace('AGT','S').replace('AGC','S').replace('TGT','0').replace('TGC','0').replace('ACT','1').replace('ACC','1').replace('ACA','1').replace('ACG','1').replace('GCT','2').replace('GCC','2').replace('GCA','2').replace('GCG','2').replace('GGT','3').replace('GGC','3').replace('GGA','3').replace('GGG','3').replace('0','C').replace('1','T').replace('2','A').replace('3','G').replace(' ','')
    

def format_fasta(seq,width=62):
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
    
    
def orfs_to_seq(orfs_list,seq_list,seq_name,out='d'):
    if not len(orfs_list) == len(seq_list):
        print ("Error")
        return ''

    ind=0
    result=[]
    for i in range(len(orfs_list)):
        
        #pair is a list [current_start_index,current_stop_index,this_frame,this_start_codon,this_stop_codon]
        pair=orfs_list[i]
        ind+=1
        #print(pair[0],pair[1])
        ostart=pair[0]
        oend=pair[1]
        frame=pair[2]+1
        startcodon=pair[3]
        stopcodon=pair[4]
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
        
        thisorfid=seq_name+"_ORF."+str(ind)+' ['+str(ostart)+'-'+str(oend)+']('+strand+') type:'+otype+' length:'+str(olen)+' frame:'+str(frame)+' start:'+startcodon+' stop:'+stopcodon
        thisseq=seq_list[i]
        if out=='p':
            #convert to prot
            thisseq=translate_dna(thisseq)
        if out == 'r':
            #conver to RNA
            thisseq=transcribe_dna(thisseq)
        result.append('>'+thisorfid+'\n'+format_fasta(thisseq))
        #print('>'+thisorfid+'\n'+seq_list[i])
    return '\n'.join(result)

#compile results in bed12
def orfs_to_bed12(orfs_list,seq_name,seqlen):
    #print(orfs_list)
    result=[]
    ind=0
    for pair in orfs_list:
        ind+=1
        #print(pair)
        #pair is a list [current_start_index,current_stop_index,this_frame,this_start_codon,this_stop_codon]
        ostart=pair[0]
        oend=pair[1]
        frame=pair[2]+1
        startcodon=pair[3]
        stopcodon=pair[4]
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
        oid= 'ID='+thisorfid+';ORF_type='+otype+';ORF_len='+str(olen)+';ORF_frame='+str(frame)+';Start:'+startcodon+';Stop:'+stopcodon
        thisorf=seq_name+'\t'+str(0)+'\t'+str(seqlen-1)+'\t'+oid+'\t'+strand+'\t'+str(ostart)+'\t'+str(oend)+'\t'+str(0)+'\t'+str(seqlen-1)+'\t'+str(0)
        #print(thisorf)
        result.append(thisorf)
    #return as string    
    return '\n'.join(result)
        
            

#TODO: implement nested,partial3 options
def get_orfs(seq,
             seqname,
             minlen,
             maxlen,
             starts=['ATG'],
             stops=['TAA','TAG','TGA'],
             nested=False, 
             partial3=False, 
             partial5=False, 
             rev_com=False):
    """

    Parameters
    ----------
    seq : TYPE
        DESCRIPTION.
    seqname : TYPE
        DESCRIPTION.
    minlen : TYPE
        DESCRIPTION.
    starts : TYPE, optional
        DESCRIPTION. The default is ['ATG'].
    stops : TYPE, optional
        DESCRIPTION. The default is ['TAA','TAG','TGA'].
    report_incomplete : TYPE, optional
        DESCRIPTION. The default is True.
    rev_com : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    return ORF positions, sequence
    allorfs : TYPE
        DESCRIPTION.
    allorfs_seq : TYPE
        DESCRIPTION.

    """
    
    cdef int seq_len=len(seq)   
    #get start and stop positions
    cdef list start_positions=[]
    cdef list stop_positions=[]
    
    #re is extremely fast
    for c in starts:
        start_positions.extend([m.start() for m in re.finditer(c,seq)])
    for c in stops:
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
    
    #if an incomplete ORF is found in a frame all upstream start codons in that frame should be ignored
    cdef int incomplete_found[3]
    incomplete_found[:]=[-1,-1,-1]
    
    for i in range(len(start_positions)):
        current_start_index=start_positions[i]
        #print('examinig',current_start_index)
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
                current_stop_found=True
                current_length=current_stop_index-current_start_index #length excluding stop codon
                #print('CL',current_length,'SI now',current_start_index)
                if current_length >= minlen:
                    #start and stop codons
                    this_stop_codon=seq[current_stop_index:current_stop_index+3]
                    this_start_codon=seq[current_start_index:current_start_index+3]
                    #Using 0 based coordinate, slice [current_start_index,current_stop_index] will yeild the ORF without the stop codon
                    complete_orfs.append([current_start_index,current_stop_index,this_frame,this_start_codon,this_stop_codon]) 
                    current_orf_seq = seq[current_start_index:current_stop_index] #current seq
                    complete_orfs_seq.append(current_orf_seq)
                    break
                else:
                    #break and forget about the ORF with smaller length
                    break
            
            
        #failed to find a stop after searching list of stops
        #if no stops left: ORF has start codon but no downstream in-frame stop codon
        if not current_stop_found:
            #print('No stops for',current_start_index)
            #this means a start codon upstream already codes for an incomplete ORF
            if incomplete_found[this_frame] < 0:
                #print('Adding',current_start_index)
                incomplete_found[this_frame]=1
                #ORFs without a stop
                if partial5:
                    current_orf_seq=seq[current_start_index:]
                    #if seq isnt multiple of 3
                    current_length=len(current_orf_seq)
                    endind=current_length-(current_length%3)
                    current_orf_seq=current_orf_seq[0:endind]
                    current_length=len(current_orf_seq)
                    #print('Addeding IC',current_start_index,' Len of IC is',current_length,'seqlen',len(current_orf_seq),current_orf_seq)
                    if current_length >= minlen:
                        this_stop_codon="NA"
                        this_start_codon=seq[current_start_index:current_start_index+3]
                        incomplete_orfs.append([current_start_index,-1,this_frame,this_start_codon,this_stop_codon])
                        incomplete_orfs_seq.append(current_orf_seq)
                        #print('Addeding IC',current_start_index)
            
                
    
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
    #return allorfs
        
        
def transform_to_sense(orfs_list,total_len):
    total_len-=1
    for orf in orfs_list:
        orf[0]=total_len-orf[0]
        if orf[1]>=0:
            orf[1]=total_len-orf[1]
        else:
            orf[1]=-9
    return orfs_list    
    


