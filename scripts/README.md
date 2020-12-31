## Helper scripts for orfipy
These scripts are helper scripts, but are not included in orfipy package.


## Running benchmarks

### Download data

- Download required data: run `bash prepare_data.sh`. This will download the benchmark data in ../testdata/bmdata directory


### Mean runtimes

The run_benchmarks.sh script runs the tools and compares the mean runtimes over n runs. 
*The runtimes are compared for each tool to report the ORFs in both nucleotide and peptide fasta format.*

To run the comparison reported in the paper:

1. Install the python package pyrpipe version 0.0.4 via `conda install -c bioconda pyrpipe` or `pip install pyrpipe`. The runtimes are recorded via pyrpipe.
2. Run the comparison:

```
 bash run_benchmarks.sh infasta.fa outdir 5 300
 
```

For example to run with human data: `bash run_benchmarks.sh ../testdata/bmdata/human.fa human_out 3 99` will run the tools 3 times each and find ORFs of length 99 and higher.
The commands executed for tools are in `run_getorf.sh`  `run_orfipy.sh` and  `run_orfm.sh` files.

### Run basic benchmark without pyrpipe
Run the basic_benchmark.sh to measure the runtime for each tool using `time`. For example:

```
bash basic_benchmark.sh ../testdata/bmdata/human.fa 99 human_basic_out
```
First argument is input fasta file, Second is the min orf length, Third is the out directory

### Reproduce figures in the paper
Use the file `makeplot.ipynb` in jupyter notebook to reproduce the figures in the paper.



### Scripts description

* `compare_fasta_files.py`: compare multiple fasta files by sequence. Example: `python compare_fasta_files.py file1.fa file2.fa file3.fa`
* `parse_NCBI_tables.py`: a python script to parse NCBI trans tables into python dict/json object. The input is in `translationtables.txt` and the output produces is `translation_tables.json`. These tables are included in orfipy.
* `run_benchmarks.sh`: runs orfipy, getorf and orfm via `pyrpipe` and produces avg runtime for each program. This script calls `benchmark.py` script which executes the tools to search regions between STOP codons and produce nucleotide and peptide fasta files. OrfM only supports this mode. Additionally, it runs getorf to search regions between start and stop codons and equivalent orfipy command. This script takes 4 arguments: 

    1. input fasta files to run the tools on
    2. output directory
    3. Num trials to run each tool
    4. Min size of ORF to find (must be multiple of 3 or OrfM fails to run). 

    Example: To run benchmarks `bash run_benchmarks.sh infasta.fa outdir 5 300`

    **Note**: `run_benchmarks.sh` requires [`pyrpipe` package](https://github.com/urmi-21/pyrpipe/) which could be installed via pip or bioconda.
