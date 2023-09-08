#!/bin/bash
# author: Bercin Cenik
# this script generates batch scripts for QC, alignment and track generation of paired-end ChIP-seq data. 
# spike-in not implemented yet.

# Initialize arguments with default values
runflag=0

# Capture user arguments
while getopts w:f:g:s:t:r: option
do
  case "${option}" in
    w) workDir=${OPTARG};;
    f) origFastqDir=${OPTARG};;
    g) genome=${OPTARG};;
    s) scientistID=${OPTARG};;
    t) tangoNumber=${OPTARG};;
    r) runflag=${OPTARG};;
    *) echo "Usage: $0 -w <workDir> -f <fastqDir> -g <genome> -s <scientistID> -t <tangoNumber> [-r <runflag: 0 or 1>]"; exit 1;;
  esac
done

# Check if all mandatory arguments are set
if [ -z "${workDir}" ] || [ -z "${origFastqDir}" ] || [ -z "${genome}" ] || [ -z "${scientistID}" ] || [ -z "${tangoNumber}" ]; then
    echo "Error: One or more mandatory options are not set!"
    echo "Usage: $0 -w <workDir> -f <fastqDir> -g <genome> -s <scientistID> -t <tangoNumber> [-r <runflag: 0 or 1>]"
    exit 1
fi


# Derived directories
fastqDir="${workDir}/fastq.clean"
fastqcDir="${workDir}/fastqc"
bamDir="${workDir}/bam"
trackDir="${workDir}/tracks"
scriptDir="${workDir}/scripts"
metadataDir="${workDir}/metadata"

# Create necessary directories
mkdir -p "${fastqDir}" "${fastqcDir}" "${bamDir}" "${trackDir}" "${scriptDir}" "${metadataDir}"

# Generate sampleNames.txt from the fastqDir
ls ${origFastqDir}/*_R1_001.fastq.gz | sed -E 's/_R1_001.fastq.gz//g' | xargs -n 1 basename > "${workDir}/sampleNames.orig.txt"
cat "${workDir}/sampleNames.orig.txt" | sed -E 's/(_S[0-9]+)?//g' > "${workDir}/sampleNames.txt"
paste "${workDir}/sampleNames.orig.txt" "${workDir}/sampleNames.txt" | while IFS=$'\t' read -r origSample sanitizedSample; do
    # Use the original sample name to construct the filename for copying
    echo "Copying sample ${origSample}..."
    cp "${origFastqDir}/${origSample}_R1_001.fastq.gz" "${fastqDir}/${sanitizedSample}_R1.fastq.gz"
    cp "${origFastqDir}/${origSample}_R2_001.fastq.gz" "${fastqDir}/${sanitizedSample}_R2.fastq.gz"
done


# Rename the files and move them to the fastq.clean directory
while read -r sample; do
    # Construct the original sample names from the sanitized one by reversing the sed operation for both R1 and R2
    origSampleName_R1=$(echo ${sample} | sed -E 's/(.*)/\1_S1_R1_001.fastq.gz/')
    origSampleName_R2=$(echo ${sample} | sed -E 's/(.*)/\1_S1_R2_001.fastq.gz/')

    # New sample names for renaming
    newSampleName_R1=$(echo ${sample} | sed -E 's/(.*)/\1_R1.fastq.gz/')
    newSampleName_R2=$(echo ${sample} | sed -E 's/(.*)/\1_R2.fastq.gz/')

    cp "${origFastqDir}/${origSampleName_R1}" "${fastqDir}/${newSampleName_R1}"
    cp "${origFastqDir}/${origSampleName_R2}" "${fastqDir}/${newSampleName_R2}"
done < "${workDir}/sampleNames.txt"

# Process each sample
while read -r sample; do
    outfile="${scriptDir}/${sample}.PE.ChIP.sh"

    cat > "${outfile}" << EOF
#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH --mem=100000
#SBATCH --chdir=${workDir}
#SBATCH -o "${metadataDir}/${sample}.o%j"
#SBATCH --error="${metadataDir}/${sample}.e%j"
#SBATCH --job-name=${tangoNumber}
#SBATCH --nodes=1
#SBATCH -n 24

# load necessary modules
module purge
module load bwa/0.7.12
module load samtools/1.6
module load pigz/2.4
module load R/3.2.2
module load picard/1.131
module load java/jdk1.8.0_25
module load bowtie2/2.2.6
module load perl/5.16
module load deeptools
module load python/anaconda

# run fastq screen for both R1 and R2
perl /projects/p20742//tools/bin/fastq_screen_v0.11.4/fastq_screen --threads 24 --aligner bowtie2 --conf /projects/p20742//tools/bin/fastq_screen_v0.11.4/fastq_screen.allRefs.conf --outdir ${fastqcDir} ${fastqDir}/${sample}_R1.fastq.gz

perl /projects/p20742//tools/bin/fastq_screen_v0.11.4/fastq_screen --threads 24 --aligner bowtie2 --conf /projects/p20742//tools/bin/fastq_screen_v0.11.4/fastq_screen.allRefs.conf --outdir ${fastqcDir} ${fastqDir}/${sample}_R2.fastq.gz

# Trim PE poor quality sequence with TRAILING:30 MINLEN:20 (see Trimmomatic documentation)
java -jar /projects/p20742//tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 24 -phred33 ${fastqDir}/${sample}_R1.fastq.gz ${fastqDir}/${sample}_R2.fastq.gz ${fastqDir}/${sample}_R1.fastq.trimmed.gz ${fastqDir}/${sample}_R1U.fastq.trimmed.gz ${fastqDir}/${sample}_R2.fastq.trimmed.gz ${fastqDir}/${sample}_R2U.fastq.trimmed.gz TRAILING:30 MINLEN:20

# Running FastQC to assess read quality before trimming.
date
/projects/p20742//tools/bin/FastQC/fastqc -o ${fastqcDir} ${fastqDir}/${sample}_R1.fastq.gz ${fastqDir}/${sample}_R2.fastq.gz
# Running FastQC to assess read quality of unpaired reads after trimming.
date
/projects/p20742//tools/bin/FastQC/fastqc -o ${fastqcDir} ${fastqDir}/${sample}_R1U.fastq.trimmed.gz ${fastqDir}/${sample}_R2U.fastq.trimmed.gz \

# Running FastQC to assess read quality of paired reads after trimming.
date
/projects/p20742//tools/bin/FastQC/fastqc -o ${fastqcDir} ${fastqDir}/${sample}_R1.fastq.trimmed.gz ${fastqDir}/${sample}_R2.fastq.trimmed.gz \
date

# genome is $genome

# BWA Alignment
bwa mem -M -t 24 -R "@RG\tID:${sample}\tSM:${sample}\tPU:nextseq\tCN:NUSeq\tPL:ILLUMINA" "/projects/p20742//anno/bwa_indexes/${genome}.fa" "${fastqDir}/${sample}_R1.fastq.trimmed.gz" "${fastqDir}/${sample}_R2.fastq.trimmed.gz" > "${bamDir}/${sample}.sam"

# Convert to BAM and mark duplicates with samblaster
/projects/p20742//tools/bin/samblaster/samblaster -i ${bamDir}/$sample.sam | samtools view -b -h -@ 24 -o ${bamDir}/$sample.unsorted.bam -

# Sort BAM file
samtools sort -@ 24 ${bamDir}/$sample.unsorted.bam -o ${bamDir}/$sample.sorted.bam

# Index BAM file
samtools index -@ 24 ${bamDir}/$sample.sorted.bam

# Mark duplicates with Picard
java -Xmx2g -jar /software/picard/1.131/picard-tools-1.131/picard.jar MarkDuplicates INPUT=${bamDir}/$sample.sorted.bam OUTPUT=${bamDir}/$sample.mdup.bam METRICS_FILE=${bamDir}/$sample.mdup_metrics.txt

# Index the marked BAM file
samtools index -@ 24 ${bamDir}/$sample.mdup.bam

# Rename the final BAM file and remove intermediates
mv ${bamDir}/$sample.mdup.bam ${bamDir}/$sample.bam
mv ${bamDir}/$sample.mdup.bam.bai ${bamDir}/$sample.bam.bai
rm ${bamDir}/$sample.sam ${bamDir}/$sample.unsorted.bam ${bamDir}/$sample.sorted.bam

# Make ChIPseq tracks.
bamCoverage --bam ${bamDir}/$sample.bam --outFileName ${trackDir}/$sample.bw --outFileFormat bigwig --extendReads --binSize 1 --scaleFactor 1 --numberOfProcessors 10

# Check if bwlist file exists, and if not, create it.
if  ! [-f "${trackDir}/bwlist.txt"];
then
echo "${trackDir}/$sample.bw" | cat > ${trackDir}/bwlist.txt 
else
echo "${trackDir}/$sample.bw" | cat >> ${trackDir}/bwlist.txt 
fi

# Make header files for tracks.
echo "track type=bigWig name=\"${sample}.bw\" description=\"${sample}.rpm\" graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/ash-tracks/TANGO/${scientistID}/${scientistID}.TANGO-${tangoNumber}/${sample}.bw" > ${trackDir}/${sample}.bw.header.txt
date

EOF

    # Check runflag
    if [ "$runflag" -eq 1 ]; then
        jobID=$(sbatch "${outfile}" | awk '{print $NF}')
        echo "Script ${outfile} has been submitted and the job ID is ${jobID}."
    else
        echo "Script ${outfile} has been generated!"
    fi

done < "${workDir}/sampleNames.txt"
