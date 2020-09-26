#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 18:57:38 2020

@author: usingh
"""
import argparse
import sys
import os
import orfipy.translation_tables
import orfipy.findorfs

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


def print_tables():
    """
    Print translation tables and exit
    """
    tab_dict=orfipy.translation_tables.translation_tables_dict
    print ('Translation tables compiled from:',tab_dict['source'])    
    print ('Table#','Name','Start','Stop')
    for k in tab_dict:
        if 'source' in k:
            continue
        print (k,
               tab_dict[k]['name'],
               '['+','.join(tab_dict[k]['start'])+']',
               '['+','.join(tab_dict[k]['stop'])+']')
    

def main():
    parser = argparse.ArgumentParser(description='orfipy: extract Open Reading Frames',
    usage="""
    orfipy [<options>] <infasta>
    If no output type, i.e. dna, rna, pep, bed or bed12, is specified, default output is bed format to stdout.
    """)
    parser.add_argument("--procs", help="Num processes\nDefault:mp.cpu_count()")
    parser.add_argument("--single-mode", help="Run in single mode i.e. no parallel processing (SLOWER). If supplied with procs, this is ignored. \nDefault: False", dest='single', action='store_true',default=False)

    parser.add_argument("--table", help="<int> The translation table to use\nDefault: 1\nLits of tables at https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes",default="1")
    parser.add_argument("--start", help="Comma-separated list of start-codons. This will override start codons described in translation table\nDefault: Derived from the translation table selected")
    parser.add_argument("--stop", help="Comma-separated list of stop codons. This will override stop codons described in translation table\nDefault: Derived from the translation table selected")
    
    #output options    
    parser.add_argument('--outdir', help='Path to outdir\ndefault: orfipy_<infasta>_out',action="store")
    parser.add_argument("--bed12", help="bed12 out file\nDefault: None")
    parser.add_argument("--bed", help="bed out file\nDefault: None")
    parser.add_argument("--dna", help="fasta (DNA) out file\nDefault: None")
    parser.add_argument("--rna", help="fasta (RNA) out file\nDefault: None")
    parser.add_argument("--pep", help="fasta (peptide) out file\nDefault: None")
    parser.add_argument("--min", help="Minimum length of ORF, excluding stop codon (nucleotide)\nDefault: 30",default=30)
    parser.add_argument("--max", help="Maximum length of ORF, excluding stop codon (nucleotide)\nDefault: inf",default=100000000000000)
    parser.add_argument("--strand", help="Strands to find ORFs [(f)orward,(r)everse,(b)oth]\nDefault: b",default='b',choices=['f', 'r', 'b'])
    parser.add_argument("--nested", help="Output nested and overlapping ORFs in the same frame \nDefault: False",default=False,dest='nested', action='store_true')
    parser.add_argument("--partial-3", help="Output ORFs lacking a start codon\nDefault: False",default=False,dest='partial3', action='store_true')
    parser.add_argument("--partial-5", help="Output ORFs lacking a stop codon\nDefault: False",default=False,dest='partial5', action='store_true')
    parser.add_argument("--longest", help="Output only longest ORFs per sequence\nDefault: False",default=False,dest='longest', action='store_true')
    parser.add_argument("--byframe", help="Write ORFs output by frame\nDefault: False",default=False,dest='byframe', action='store_true')
    
    
    
    parser.add_argument("--chunk-size", help="""Max chunk size in MB. This is useful for limiting memory usage when processing large fasta files using multiple processes
                        The files are processed in chunks if file size is greater than chunk-size. 
                        By default orfipy computes the chunk size based on available memory and cpu cores.
                        Providing a smaller chunk-size will lower the memory usage but, actual memory used by orfipy can be more than the chunk-size. 
                        Providing a very high chunk-size can lead to memory problems. 
                        It is best to let orfipy decide on the chunk-size.
                        \nDefault: estimated by orfipy based on system memory and cpu""")
                        
    parser.add_argument("--show-tables", help="Print translation tables and exit.\nDefault: False",default=False,dest='showtab', action='store_true')
    
    
    parser.add_argument('infile', help='Fasta file with input DNA sequences',action="store",nargs="?")
    
    
    args = parser.parse_args()
    
    
    
    if args.showtab:
        print_tables()
        sys.exit(0)
        
    infile=args.infile    
    if not infile:
        parser.print_help()
        sys.exit(1)
    
    #parse args
    minlen=args.min
    if minlen:
        minlen=int(minlen)
    else:
        minlen=30
    #maxlen
    maxlen=args.max
    if maxlen:
        maxlen=int(maxlen)
    else:
        maxlen=1000000000000
    
    
    single=args.single
    procs=args.procs
    if procs:
        procs=int(procs)
        single=False
        
    #init translation table
    tablenum=args.table
    table=orfipy.translation_tables.translation_tables_dict[tablenum]
    starts=table['start']
    stops=table['stop']
    #if start and stop are provided, override
    if args.start:
        starts=args.start
        starts=starts.split(',')
    if args.stop:
        stops=args.stop
        stops=stops.split(',')
    #valid start, stop codons
    if not validate_codons(starts,stops):
        print('Please check start/stop codon list again')
        sys.exit(1)
    print('Using translation table:',table['name'],'start:',starts,'stop:',stops)
        
    
    strand=args.strand
    #chunk-size; if not provided, it is estimated in findorfs.py
    chunk_size=args.chunk_size
    
    bed12=args.bed12
    bed=args.bed
    dna=args.dna
    rna=args.rna
    pep=args.pep
    
    outdir=args.outdir
    if not outdir:
        outdir="orfipy_"+os.path.basename(infile)+'_out'
        #print("Temp dir is {}".format(tmpdir),file=sys.stderr)
    
    
    #print(args)
    #call main program    
    orfipy.findorfs.main(infile,minlen,maxlen,procs,single,chunk_size,strand,starts,stops,args.nested,args.partial3,args.partial5,args.longest,args.byframe,bed12,bed,dna,rna,pep,outdir)
    
    
    
if __name__ == '__main__':
    main()
