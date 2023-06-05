# Draft of ONT bam merge scripts

Created by: Bercin Cenik
Date: June 5, 2023 10:51 AM

Overview: this is a list of all the steps I perform for Oxford Nanopore cDNA bam file merging. It also contains instructions for indexing the files with `samtools`

### Merging of alignments:

Typically, PromethION runs generate 5,000-15,000 minibam files which will need to be merged. Since merging them all at once can be computationally too intensive even on Quest, I do this in 4 successive steps (R1-R2-R3-final) that iteratively merge the files in groups of 10 using samtools.

The scripts for creating the merge commands are:

`create.merge.scripts.sh`

`create.R2.scripts.sh`

`create.R3.scripts.sh`

First, set your file names and directories:

```bash
directory=/projects/b1042/Shilatifard/PromethION/CHACHAx/minibam
subdirectory=bam_pass/
file1='file1'
file2='file2'
file3='file3'
file4='file4'
```

Include this in all shell scripts that generate the merge scripts.

Set up your directory as follows

CHACHAx

fastq

seq.summary

qc.reports

tars

bam

minibam

nameofsample

bam_pass

minibam files

nameofsample2

bam_pass

minibam files

etc.

create.merge.scripts.sh

```bash
## make the file path lists
## this generates a text file with the full paths to all minibam files

for file in "file1" "file2" "file3" "file4"
do
ls ${directory}/${!file}/${subdirectory}/*.bam > ${!file}.minibamlist.txt
done

## create the first merge script to generate the R1 bams:
for item in "file1" "file2" "file3" "file4"
do
    touch ${!item}.scripts.sh
    echo "#!/bin/sh" > ${!item}.scripts.sh
    echo "# Bercin Cenik" >> ${!item}.scripts.sh
    echo "# Bercin Cenik" >> ${!item}.scripts.sh
    echo "#SBATCH -A b1025" >> ${!item}.scripts.sh
    echo "#SBATCH -p buyin" >> ${!item}.scripts.sh
    echo "#SBATCH -t 24:00:00" >> ${!item}.scripts.sh
    echo "#SBATCH -N 1" >> ${!item}.scripts.sh
    echo "#SBATCH --mem=64G" >> ${!item}.scripts.sh
    echo "#SBATCH --ntasks-per-node=10" >> ${!item}.scripts.sh
    echo "#SBATCH -o %x.o%j" >> ${!item}.scripts.sh
    echo "#SBATCH -e %x.e%j" >> ${!item}.scripts.sh
    echo " " >> ${!item}.scripts.sh   
    echo " " >> ${!item}.scripts.sh
    echo "module load samtools" >> ${!item}.scripts.sh
    echo " " >> ${!item}.scripts.sh
    sort -V ${!item}.minibamlist.txt > ${!item}.minibamlist.ordered.txt
    cat ${!item}.minibamlist.ordered.txt | xargs -n 10 > ${!item}.minibamlist.ordered.grouped.txt
    line_count=$(wc -l < ${!item}.minibamlist.ordered.grouped.txt)
    echo "Number of lines in ${!item}.minibamlist.ordered.grouped.txt: ${line_count}"
    seq -f "samtools merge -@ 8 ./merge.R1/merged_R1_%04g.bam" 1 ${line_count} > "${!item}.headers.txt"
    paste -d ' ' ${!item}.headers.txt ${!item}.minibamlist.ordered.grouped.txt >> ${!item}.scripts.sh
    ##clean up files
    #rm ${!item}.minibamlist.txt
    rm ${!item}.minibamlist.ordered.txt
    rm ${!item}.minibamlist.ordered.grouped.txt
    rm ${!item}.headers.txt
    ##make directories and subdirectories
    mkdir merge.${!item}
    mkdir merge.${!item}/merge.R1
    mkdir merge.${!item}/merge.R2
    mkdir merge.${!item}/merge.R3
    mv ${!item}.scripts.sh merge.${!item}
    echo "done with ${!item}"
done
```

This will output nameofsample.scripts.sh and move it to a directory called merge.nameofsample (different from nameofsample where the original minibams are).

Run this script on Quest, wait for it to finish. Meanwhile, create the R2 merging scripts:

create.R2.scripts.sh

```bash
# make the grouped R1 file lists
for item in "file1" "file2" "file3" "file4"
do
directory=/projects/b1042/Shilatifard/PromethION/CHACHA4/minibam/merge.${!item}
cat ${directory}/${!item}.scripts.sh | tail -n +14 | cut -f 5 -d " " > ${directory}/${!item}.R1list.txt
sort -V ${directory}/${!item}.R1list.txt > ${directory}/${!item}.R1list.ordered.txt
cat ${directory}/${!item}.R1list.ordered.txt | xargs -n 10 > ${directory}/${!item}.R1list.ordered.grouped.txt
done

# make the merge scripts
for item in "file1" "file2" "file3" "file4"
do
    directory=/projects/b1042/Shilatifard/PromethION/CHACHA4/minibam/merge.${!item}
    cd ${directory}
    touch ${!item}.mergeR2scripts.sh
    echo "#!/bin/sh" >> ${!item}.mergeR2scripts.sh
    echo "# Bercin Cenik" >> ${!item}.mergeR2scripts.sh
    echo "#SBATCH -A b1025" >> ${!item}.mergeR2scripts.sh
    echo "#SBATCH -p buyin" >> ${!item}.mergeR2scripts.sh
    echo "#SBATCH -t 24:00:00" >> ${!item}.mergeR2scripts.sh
    echo "#SBATCH -N 1" >> ${!item}.mergeR2scripts.sh
    echo "#SBATCH --mem=64G" >> ${!item}.mergeR2scripts.sh
    echo "#SBATCH --ntasks-per-node=10" >> ${!item}.mergeR2scripts.sh
    echo "#SBATCH -o %x.o%j" >> ${!item}.mergeR2scripts.sh
    echo "#SBATCH -e %x.e%j" >> ${!item}.mergeR2scripts.sh
    echo " " >> ${!item}.mergeR2scripts.sh
    echo "module load samtools" >> ${!item}.mergeR2scripts.sh
    echo " " >> ${!item}.mergeR2scripts.sh
    line_count=$(wc -l < ${directory}/${!item}.R1list.ordered.grouped.txt) 
    echo "Number of lines in ${!item}.R1list.ordered.grouped.txt: ${line_count}"
    seq -f "samtools merge -@ 8 ./merge.R2/merged_R2_%03g.bam" 1 ${line_count} > "${!item}.mergeR2.headers.txt"
    paste -d ' ' ${!item}.mergeR2.headers.txt ${directory}/${!item}.R1list.ordered.grouped.txt >> ${!item}.mergeR2scripts.sh
done
```

This outputs nameofsample.mergeR2scripts.sh into the merge.nameofsample directory, where you can submit it to the scheduler once *.scripts.sh is done running.

The R3 script generation is much of the same:

create.R3.scripts.sh

```bash
# make the grouped R2 file lists in each subdirectory
for item in "file1" "file2" "file3" "file4"
do
directory=/projects/b1042/Shilatifard/PromethION/CHACHA4/minibam/merge.${!item}
cat ${directory}/${!item}.mergeR2scripts.sh | tail -n +14 | cut -f 5 -d " " > ${directory}/${!item}.R2list.txt
sort -V ${directory}/${!item}.R2list.txt > ${directory}/${!item}.R2list.ordered.txt
cat ${directory}/${!item}.R2list.ordered.txt | xargs -n 10 > ${directory}/${!item}.R2list.ordered.grouped.txt
done

# make the merge scripts
for item in "file1" "file2" "file3" "file4"
do
    directory=/projects/b1042/Shilatifard/PromethION/CHACHA4/minibam/merge.${!item}
    touch ${directory}/${!item}.mergeR3scripts.sh
    echo "#!/bin/sh" >> ${directory}/${!item}.mergeR3scripts.sh
    echo "# Bercin Cenik" >> ${directory}/${!item}.mergeR3scripts.sh
    echo "#SBATCH -A b1025" >> ${directory}/${!item}.mergeR3scripts.sh
    echo "#SBATCH -p buyin" >> ${directory}/${!item}.mergeR3scripts.sh
    echo "#SBATCH -t 24:00:00" >> ${directory}/${!item}.mergeR3scripts.sh
    echo "#SBATCH -N 1" >> ${directory}/${!item}.mergeR3scripts.sh
    echo "#SBATCH --mem=64G" >> ${directory}/${!item}.mergeR3scripts.sh
    echo "#SBATCH --ntasks-per-node=10" >> ${directory}/${!item}.mergeR3scripts.sh
    echo "#SBATCH -o %x.o%j" >> ${directory}/${!item}.mergeR3scripts.sh
    echo "#SBATCH -e %x.e%j" >> ${directory}/${!item}.mergeR3scripts.sh
    echo " " >> ${directory}/${!item}.mergeR3scripts.sh
    echo "module load samtools" >> ${directory}/${!item}.mergeR3scripts.sh
    echo " " >> ${directory}/${!item}.mergeR3scripts.sh
    line_count=$(wc -l < ${directory}/${!item}.R2list.ordered.grouped.txt)
    echo "Number of lines is ${line_count}"
    seq -f "samtools merge -@ 8 ./merge.R3/merged_R3_%03g.bam" 1 ${line_count} > "${directory}/${!item}.mergeR3.headers.txt"
    paste -d ' ' ${directory}/${!item}.mergeR3.headers.txt ${directory}/${!item}.R2list.ordered.grouped.txt >> ${directory}/${!item}.mergeR3scripts.sh
done
```

After all files are merged, transfer the final bam files to the bam directory. The intermediate files can be deleted. This can also be incorporated into the original script so it won’t have to be done manually, but I prefer not to do it in case the run crashes.

### Sort and index the bam files

Ensure that you use at least 32G of memory and 7 sort cores. Sorting is very memory-intensive, a bam file >100GB will take at least 3 hours to sort with 64GB of memory. Indexing is incrementally less memory intensive and will take about 45-60 minutes for a similarly sized file.

sample.sort.sh

```bash
#!/usr/bin/sh
# Bercin Cenik
#SBATCH -A b1025
#SBATCH -p buyin
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=10
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
module load samtools

outDir=/projects/b1042/Shilatifard/PromethION/CHACHAx/bam
filename=filename

date
echo "sort starts"
samtools sort -@ 7 -m 4G -o ${outDir}/$filename.sorted.bam ${outDir}/$filename.bam
echo "sort ends"
date
echo "indexing starts"
samtools index -@ 8 ${outDir}/$filename.sorted.bam 
echo "indexing ends"
date
```

Afterwards, run samtools flagstat:

flagstat.sh

```bash
#!/bin/sh
# Bercin Cenik
#SBATCH -A b1025
#SBATCH -p buyin
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=10
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH --job-name=flagstat

module load samtools

outDir=/projects/b1042/Shilatifard/PromethION/CHACHA4/bam
file1=file1
file2=file2
file3=file3
file4=file4

for file in $file1 $file2 $file3 $file4
do
    samtools flagstat ${outDir}/$file.sorted.bam > $file.flagstat.txt
done
```
