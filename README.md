## ePyGenetics-Tools

* This is a collection of scripts I use for molecular biology.
### sgRNA tools
1. RevComp: [biopython](https://biopython.org/) wrapper script that allows the user to input a short nucleotide sequence and outputs the reverse complement of the sequence.
2. GuideToOligo: takes a list of 20-bp sgRNA sequences and ouputs a list of cloning ready oligos for each. Also utilizes `biopython`.
3. `get_reverse_complements.ipynb`: this jupyter notebook has a few python scripts for outputting reverse complemented oligo tables.
4. retrieving sgRNAs from library: this is a markdown file with my instructions for retrieving a set of sgRNAs from a library and formatting them into order-ready oligos.
5. ChopDatUp: breaks sequences larger than 2kb into smaller chunks for submitting to [CRISPick](https://portals.broadinstitute.org/gppx/crispick/public/). Run using following command: `python ChopDatUp.v5.py <path to input files> <path to output files>`.

### Oxford Nanopore Sequencing tools
1. `make.tar.sh` and `make.merged.fastq.sh`: these scripts respectively tarball ONT minibam files and merge ONT minifastq files. the instructions are contained in the `maketar_makemergedfastq.instructions.md` file.
2. `babam.md`: draft of a script collection for creating merged bam files from ONT minibam files, aka big-a** bam files (BA-BAM files). UPDATE: Deprecated. Please use the master script `babam.sh` instead, see below.

### Pipelines 
1. `make.PE.chip.scripts.sh`: generates individual shell scripts for QC, trimming, alignment, bigwig track generation and upload. Usage: `make.PE.chip.scripts.sh -w <workDir> -f <fastqDir> -g <genome> -s <scientistID> -t <tangoNumber> [-r <runflag: 0 or 1>]`. Note that this assumes the `_SXX_R1_001.fastq.gz` naming convention for fastq files. Modify accordingly if not the case. Spike-in version is being tested.
2. `babam.sh`: generates sequential merge scripts for minibam files generated during ONT sequencing. Usage: `bash babam.sh $0 -m <text file with absolute paths to minibam files> -w <absolute path to working directory> -r <1|0 user run flag>`. In progress!
