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
    orfipy [<options>] <infasta> <outdir>
     Specify at least one output type i.e. dna, rna, pep, bed or bed12. 
     If not specified output is bed format to stdout.
    """)
    parser.add_argument("--procs", help="Num processes\nDefault:mp.cpu_count()")
    parser.add_argument("--single-mode", help="Run in single mode i.e. no parallel processing (SLOWER). If supplied with procs, it is ignored. \nDefault: False", dest='single', action='store_true',default=False)

    
    parser.add_argument("--starts", help="Comma-separated list of start-codons\nDefault: ATG")
    parser.add_argument("--stops", help="Comma-separated list of stop codons\nDefault: TAA,TGA,TAG")
    
    #output options    
    parser.add_argument("--bed12", help="bed12 out file\nDefault: None")
    parser.add_argument("--bed", help="bed out file\nDefault: None")
    parser.add_argument("--dna", help="fasta (DNA) out file\nDefault: None")
    parser.add_argument("--rna", help="fasta (RNA) out file\nDefault: None")
    parser.add_argument("--pep", help="fasta (peptide) out file\nDefault: None")
    
    parser.add_argument("--min", help="Minimum length of ORF\ndefault: 30'",default=30)
    parser.add_argument("--strand", help="Strands to find ORFs [(f)orward,(r)everse,(b)oth]\ndefault: b",default='b',choices=['f', 'r', 'b'])
    
    parser.add_argument("--chunk-size", help="Max chunk size in MB. This is useful for limiting memory usage for large fasta files. The files are processed in chunks if file size is greater than chunk. NOTE: Having smaller chunk size lowers memory usage but chunk, actual memory used by orfipy can be double the chunk size. For python < 3.8 this can not be more that 2000\nDefault: based on system memory and py version")
    
    parser.add_argument('--outdir', help='Path to outdir',action="store")
    parser.add_argument('infile', help='Fasta containing sequences',action="store")
    
    
    args = parser.parse_args()
    
    infile=args.infile
    
    #parse args
    minlen=args.min
    if minlen:
        minlen=int(minlen)
    else:
        minlen=30
    single=args.single
    procs=args.procs
    if procs:
        procs=int(procs)
        single=False
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
    
    strand=args.strand
    chunk_size=args.chunk_size
    
    bed12=args.bed12
    bed=args.bed
    dna=args.dna
    rna=args.rna
    pep=args.pep
    
    outdir=args.outdir
    
    #call main program    
    orfipy.findorfs.main(infile,minlen,procs,single,chunk_size,strand,starts,stops,bed12,bed,dna,rna,pep,outdir)

      
    
    
    
if __name__ == '__main__':
    main()