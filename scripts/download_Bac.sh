#!/bin/bash

#inspired from https://www.biostars.org/p/61081/

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt

head -5000 assembly_summary.txt | awk -F '\t' '{if($12=="Complete Genome") print $20}' > assembly_summary_complete_genomes.txt
mkdir tmpbac
for next in $(cat assembly_summary_complete_genomes.txt); do wget -P tmpbac "$next"/*genomic.fna.gz; done
gunzip tmpbac/*.gz

cat tmpbac/*.fna > bac5000.fasta

