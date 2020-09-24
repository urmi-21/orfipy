#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: usingh
"""
import orfipy
import orfipy.findorfs as fo
import orfipy.orfipy_core as oc

seq='ATGTTTATGAAATAGAACTAAATGCCCATG'
orfs_pos=[[0, 12, 0, 'ATG', 'TAG', 'complete', 12], 
         [6, 12, 0, 'ATG', 'TAG', 'nested', 6], 
         [15, 18, 0, 'AAC', 'TAA', '3-prime-partial', 3], 
         [21, -1, 0, 'ATG', 'NA', '5-prime-partial', 9], 
         [-3, 7, 1, 'NA', 'TGA', '3-prime-partial', 6]]
orfs_seq=['ATGTTTATGAAA', 'ATGAAA', 'AAC', 'ATGCCCATG', 'TGTTTA']
         
def test_orf_search():
    result=oc.get_orfs(seq,"seq1",0,10000,rev_com=False,nested=True, partial3=True, partial5=True)
    print(result)
    assert result[0]==orfs_pos
    assert result[1]==orfs_seq

def test_transcribe():
    res=oc.transcribe_dna('ATGCTTGA')
    print(res)
    assert res=='AUGCUUGA'
    
def test_translate():
    res=oc.translate_dna('ATGCCCCTTGAGTAG')
    print(res)
    assert res=='MPLE.'