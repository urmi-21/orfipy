#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 18:57:38 2020

@author: usingh
"""
import argparse
import orfipy.findorfs
import sys

def validate_codons(starts,stops):
    validalphabets=['A','C','T','G']
    if not starts:
        starts=[]
    if not stops:
        stops=[]
    #check lengths
    for c in starts:
        if not len(c)==3:
            print('Invalid start codon:'+c)
            return False
        if not ((c[0] in validalphabets) and (c[1] in validalphabets) and (c[2] in validalphabets)):
            print('Invalid start codon:'+c)
            return False
    for c in stops:
        if not len(c)==3:
            print('Invalid stop codon:'+c)
            return False
    #check intersection
    for c in starts:
        if c in stops:
            print('Invalid codon in both lists:'+c)
            return False
    return True

def main():
    parser = argparse.ArgumentParser(description='orfipy: extract Open Reading Frames',
    usage="""
    orfipy [<options>] <infasta>
    """)
    parser.add_argument("--procs", help="Num processes\nDefault:mp.cpu_count()")
    parser.add_argument("--starts", help="Comma-separated list of start-codons\nDefault: ATG")
    parser.add_argument("--stops", help="Comma-separated list of stop codons\nDefault: TAA,TGA,TAG")
    parser.add_argument("--out", help="Path to out file\nDefault: None (stdout)")
    parser.add_argument("--outfmt", help="Outfile type [(b)ed,(d)na,(r)na,(p)eptide]\ndefault: a'",default='b',choices=['b', 'd', 'r','p'])
    parser.add_argument("--min", help="Minimum length of ORF\ndefault: 30'",default=30)
    parser.add_argument('infile', help='Fasta containing sequences',action="store")
    args = parser.parse_args()
    
    infile=args.infile
    
    #parse args
    minlen=args.min
    if minlen:
        minlen=int(minlen)
    procs=args.procs
    if procs:
        procs=int(procs)
    starts=args.starts
    if starts:
        starts=starts.split(',')
    else:
        starts=['ATG']
    stops=args.stops
    if stops:
        stops=stops.split(',')
    else:
        stops=['TAA','TAG','TGA']
    if not validate_codons(starts,stops):
        print('Please check start/stop codon list again')
        sys.exit(1)
    outfile=args.out
    outfmt=args.outfmt
        
    #print(args)
    #print(minlen,procs,starts,stops,outfile,outfmt)
    
    #call main program    
    orfipy.findorfs.main(infile,minlen,procs,starts,stops,outfile,outfmt)

      
    
    
    
if __name__ == '__main__':
    main()