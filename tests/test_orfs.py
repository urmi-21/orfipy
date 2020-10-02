#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: usingh
"""
import orfipy
import orfipy.findorfs as fo
import orfipy.orfipy_core as oc

seq='ATGTTTATGAAATAGAACTAAATGCCCATG'
orfs_pos1=[(0, 12, 1, 'ATG', 'TAG', 'complete', 12), 
           (1, 7, 2, 'TGT', 'TGA', '5-prime-partial', 6), 
           (2, 29, 3, 'GTT', 'NA', '5-prime-partial', 27), 
           (10, 28, 2, 'AAT', 'NA', '5-prime-partial', 18), 
           (15, 18, 1, 'AAC', 'TAA', '5-prime-partial', 3), 
           (21, 30, 1, 'ATG', 'NA', '3-prime-partial', 9)] 
orfs_seq1=['ATGTTTATGAAA', 'TGTTTA', 'GTTTATGAAATAGAACTAAATGCCCAT', 'AATAGAACTAAATGCCCA', 'AAC', 'ATGCCCATG']


def test_orf_search():
    result=oc.get_orfs(seq,"seq1",0,10000,rev_com=False, partial3=True, partial5=True,return_seqs=True)
    print(result)
    assert result[0]==orfs_pos1
    assert result[1]==orfs_seq1

