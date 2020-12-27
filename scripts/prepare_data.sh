#!/bin/bash

mkdir -p ../testdata/bmdata

wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O ../testdata/bmdata/human.fa.gz
gunzip --force ../testdata/bmdata/human.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-48/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -O ../testdata/bmdata/at.fa.gz
gunzip --force ../testdata/bmdata/at.fa.gz

cp download_Bac.sh ../testdata/bmdata/
cd ../testdata/bmdata/
bash download_Bac.sh

