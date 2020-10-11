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
import json

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
    orfipy [<options>] <infile>
    By default orfipy reports ORFs as sequences between start and stop codons. See ORF searching options to change this behaviour.
    If no output type, i.e. dna, rna, pep, bed or bed12, is specified, default output is bed format to stdout.
    """)
    parser.add_argument("--procs", help="Num processes\nDefault:mp.cpu_count()")
    parser.add_argument("--single-mode", help="Run in single mode i.e. no parallel processing (SLOWER). If supplied with procs, this is ignored. \nDefault: False", dest='single', action='store_true',default=False)

    parser.add_argument("--table", help="The codon table number to use or path to .json file with codon table.\nUse --show-tables to see available tables compiled from: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes\nDefault: 1",default="1")
    parser.add_argument("--start", help="Comma-separated list of start-codons. This will override start codons described in translation table. E.g. \"--start ATG,ATT\"\nDefault: Derived from the translation table selected")
    parser.add_argument("--stop", help="Comma-separated list of stop codons. This will override stop codons described in translation table. E.g. \"--start TAG,TTT\"\nDefault: Derived from the translation table selected")
    
    #output options    
    parser.add_argument('--outdir', help='Path to outdir\ndefault: orfipy_<infasta>_out',action="store")
    parser.add_argument("--bed12", help="bed12 out file\nDefault: None")
    parser.add_argument("--bed", help="bed out file\nDefault: None")
    parser.add_argument("--dna", help="fasta (DNA) out file\nDefault: None")
    parser.add_argument("--rna", help="fasta (RNA) out file\nDefault: None")
    parser.add_argument("--pep", help="fasta (peptide) out file\nDefault: None")
    parser.add_argument("--min", help="Minimum length of ORF, excluding stop codon (nucleotide)\nDefault: 30",default=30)
    parser.add_argument("--max", help="Maximum length of ORF, excluding stop codon (nucleotide)\nDefault: 1,000,000,000",default=1000000000)
    parser.add_argument("--strand", help="Strands to find ORFs [(f)orward,(r)everse,(b)oth]\nDefault: b",default='b',choices=['f', 'r', 'b'])
    #parser.add_argument("--nested", help="Output nested and overlapping ORFs in the same frame \nDefault: False",default=False,dest='nested', action='store_true')
    parser.add_argument("--partial-3", help="Output ORFs with a start codon but lacking an inframe stop codon. E.g. \"ATG TTT AAA\"\nDefault: False",default=False,dest='partial3', action='store_true')
    parser.add_argument("--partial-5", help="Output ORFs with an inframe stop codon lacking an inframe start codon. E.g. \"TTT AAA TAG\"\nDefault: False",default=False,dest='partial5', action='store_true')
    parser.add_argument("--between-stops", help="Output ORFs defined as regions between stop codons (regions free of stop codon). This will set --partial-3 and --partial-5 true.\nDefault: False",default=False,dest='bw_stops', action='store_true')
    parser.add_argument("--include-stop", help="Include stop codon in the results, if a stop codon exists. This output format is compatible with TransDecoder's which includes stop codon coordinates\nDefault: False",default=False,dest='include_stop', action='store_true')
    parser.add_argument("--longest", help="Output a separate BED file for longest ORFs per sequence. Requires bed option.\nDefault: False",default=False,dest='longest', action='store_true')
    parser.add_argument("--by-frame", help="Output separate BED files for ORFs by frame. Can be combined with \"--longest\" to output longest ORFs in each frame. Requires bed option.\nDefault: False",default=False,dest='byframe', action='store_true')
    
    
    
    parser.add_argument("--chunk-size", help="""Max chunk size in MB. This is useful for limiting memory usage when processing large fasta files using multiple processes
                        The files are processed in chunks if file size is greater than chunk-size. 
                        By default orfipy computes the chunk size based on available memory and cpu cores.
                        Providing a smaller chunk-size will lower the memory usage but, actual memory used by orfipy can be more than the chunk-size. 
                        Providing a very high chunk-size can lead to memory issues for larger sequences such as large chromosomes. 
                        It is best to let orfipy decide on the chunk-size.
                        \nDefault: estimated by orfipy based on system memory and cpu""")
                        
    parser.add_argument("--show-tables", help="Print translation tables and exit.\nDefault: False",default=False,dest='showtab', action='store_true')
    
    
    parser.add_argument('infile', help='The input file, in Fasta format, containing Nucletide sequences',action="store",nargs="?")
    
    
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
        maxlen=1000000000
    
    
    single=args.single
    procs=args.procs
    if procs:
        procs=int(procs)
        single=False
        
    #init translation table
    tablenum=args.table
    #check if file is specified as table
    if tablenum.isnumeric():
        table=orfipy.translation_tables.translation_tables_dict[tablenum]
        starts=table['start']
        stops=table['stop']
    else:
        if not os.path.exists(tablenum):
            print('Invalid Table Provided:',tablenum)
            sys.exit(1)
        #load table
        with open(tablenum) as f:
            table = json.load(f)
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
    print('Using translation table:',table['name'],'start:',starts,'stop:',stops,file=sys.stderr)
        
    
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
        
    #if longest of byframe is specified make sure bed format is on
    if args.longest or args.byframe:
        if not bed:
            print('Please specify the bed output option if providing --longest or --byframe')
            sys.exit(1)
    
    #print(args)
    #call main program    
    #print('VersionXXXXX')
    orfipy.findorfs.main(infile,
                         minlen,
                         maxlen,
                         procs,
                         single,
                         chunk_size,
                         strand,
                         starts,
                         stops,
                         table['table'],
                         args.include_stop,
                         args.partial3,
                         args.partial5,
                         args.bw_stops, 
                         args.longest,
                         args.byframe,
                         bed12,
                         bed,
                         dna,
                         rna,
                         pep,
                         outdir)
    
    
    
if __name__ == '__main__':
    main()
