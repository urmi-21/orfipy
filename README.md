[![Build Status](https://travis-ci.org/urmi-21/orfipy.svg?branch=master)](https://travis-ci.org/urmi-21/orfipy)

# Introduction
orfipy is a python package to extract ORFs in extremely fast and flexible manner. 


## Installation
Install latest stable version
```
pip install orfipy
```

Install development version source

```
git clone https://github.com/urmi-21/orfipy.git
cd orfipy
pip install .
```

## Examples

Extract ORF sequences and write ORF sequences in orfs.fa file

```
orfipy input.fasta --dna orfs.fa --min 10 --max 10000 --procs 4 --table 1 --outdir orfs_out
```

Use [standard codon table](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes)  but use only ATG as start codon

```
orfipy input.fasta --dna orfs.fa --start ATG
```

Extract ORF BED file
```
orfipy input.fasta --bed orfs.bed --min 50 --procs 4
or
orfipy input.fasta --min 50 --procs 4 > orfs.bed 
```

Extract ORF BED12 file

**Note**: Add `--include-stop` for orfipy output to be consistent with Transdecoder.Predict output .bed file

```
orfipy testseq.fa --min 100 --bed12 of.bed --partial-5 --partial-3 --include-stop
```

Extract peptides
```
orfipy input.fasta --pep orfs_peptides.fa --min 50 --procs 4
```

See available codon tables
```
orfipy --show-table

```

## Comparison with Getorf and OrfM

![alt text](https://raw.githubusercontent.com/urmi-21/orfipy/master/scripts/benchmarks.png)














