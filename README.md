# Introduction
orfipy is a python package to extract ORFs in extremely fast and flexible manner.

## Stability
Experimental

## Installation
Install from source

```
git clone 
pip install -e <path to orfipy>
```

## Examples

Extract ORF sequences and write ORF sequences in orfs.fa file
```
orfipy input.fasta --dna orfs.fa --min 10 --max 10000 --procs 4 --table 1
```

Extract ORF BED file
```
orfipy input.fasta --bed orfs.bed --min 50 --procs 4
or
orfipy input.fasta --min 50 --procs 4 > orfs.bed 
```

Extract peptides
```
orfipy input.fasta --pep orfs_peptides.fa --min 50 --procs 4
```
