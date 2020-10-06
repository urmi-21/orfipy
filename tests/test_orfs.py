#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: usingh
"""
import orfipy_core as oc

seq='ATGTTTATGAAATAGAACTAAATGCCCATG'
seq_rc='CATGGGCATTTAGTTCTATTTCATAAACAT'


#start_search(seq,seq_rc,seqname,minlen,maxlen,strand,starts,stops,table,include_stop,partial3,partial5,find_between_stops,out_opts)
def test_orf_search():
    result=oc.start_search(seq,
                    seq_rc,
                    'test',
                    0,
                    100000,
                    'b',
                    ['ATG'],
                    ['TAG'],
                    '1',
                    True,
                    True,
                    False,
                    False,
                    [False,False,False,False,False]       
        
        )
    
    print(result)
    print(result[0])
    print(len(result[0].split('\n')))
    assert len(result[0].split('\n'))==3
    
    