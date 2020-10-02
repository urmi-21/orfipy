### Helper scripts for orfipy
These scripts are helper scripts, but are not included in orfipy package.
Important scripts are:

* `parse_NCBI_tables.py`: a python script to parse NCBI trans tables into python dict/json object. The input is in `translationtables.txt` and the output produces is `translation_tables.json`. These tables are included in orfipy.
* `run_benchmarks.sh`: runs orfipy, getorf and orfm via `pyrpipe` and produces avg runtime for each program. This script calls `benchmark.py` script which executes the tools to search regions between STOP codons and produce nucleotide and peptide fasta files. OrfM only supports this mode. Additionally, it runs getorf to search regions between start and stop codons and equivalent orfipy command. This script takes 4 arguments: 
    1. input fasta files to run the tools on
    2. output directory
    3. Num trials to run each tool
    4. Min size of ORF to find (must be multiple of 3 or OrfM fails to run)
Example: To run benchmarks `bash run_benchmarks.sh infasta.fa outdir 5 30`
