## ePyGenetics-Tools
* This is a collection of scripts I use for molecular biology.
* Currently includes:
1. RevComp: [biopython](https://biopython.org/) wrapper script that allows the user to input a short nucleotide sequence and outputs the reverse complement of the sequence.
2. GuideToOligo: takes a list of 20-bp sgRNA sequences and ouputs a list of cloning ready oligos for each. Also utilizes `biopython`.
3. get_reverse_complements: this jupyter notebook has a few python scripts for outputting reverse complemented oligo tables.
4. retrieving sgRNAs from library: this is a markdown file with my instructions for retrieving a set of sgRNAs from a library and formatting them into order-ready oligos.
5. make.tar.sh and make.merged.fastq.sh: these scripts respectively tarball ONT minibam files and merge ONT minifastq files. the instructions are contained in the maketar_makemergedfastq.instructions.md file.
