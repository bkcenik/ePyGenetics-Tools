#!/bin/bash
# BA-BAM! A parallelized script for sequential merging of ONT minibam files
# by Bercin Cenik
# V1.0 as of 09122023
# Notes: while I optimize this script, please continue to use absolute paths.
# for the minibam files, make sure the paths are provided as: <sample name>/<ont run ID>/bam_pass
# if you want to keep the intermediate files, comment out the clean-up steps.
# Next version will incorporate sorting and indexing

# Function to print usage
print_usage() {
    echo "Usage: bash babam.sh $0 -m <text file with absolute paths to minibam files> -w <absolute path to working directory> -r <1|0 user run flag>"
}

# Process arguments
while getopts "m:w:r:" opt; do
    case $opt in
        m) bampaths="$OPTARG";;
        w) working_directory="$OPTARG";;
        r) run_flag="$OPTARG";;
        \?) print_usage
            exit 1 ;;
    esac
done

# Check if all arguments are provided
if [ -z "$bampaths" ] || [ -z "$working_directory" ] || [ -z "$run_flag" ]; then
    print_usage
    exit 1
fi

# Define log file
logfile="${working_directory}/babam.log"

# Start logging
exec &> >(tee -a "$logfile")

cd "$working_directory" || { echo "Failed to change to directory $working_directory"; exit 1; }


echo "BA-BAM! A series of scripts for merging ONT minibams.."

while read -r path; do
    item=$(echo "$path" | awk -F'/' '{print $(NF-2)}')
    sample_directory="${working_directory}/merge.${item}"
    echo "item name is $item"
    echo "sample directory is ${sample_directory}"

    # Make necessary directories
    mkdir -p "${sample_directory}/merge.R1" "${sample_directory}/merge.R2" "${sample_directory}/merge.R3" "${sample_directory}/merge.R4"
    ls "${path}"/*.bam > "$item.minibamlist.txt"

    # Create a SLURM script for merging
    filename="${item}.scripts.sh"
    echo "#! /bin/bash" > "$filename"
    echo "# Bercin Cenik" >> "$filename"
    echo "#SBATCH -A b1042" >> "$filename"
    echo "#SBATCH -p genomics" >> "$filename"
    echo "#SBATCH -t 24:00:00" >> "$filename"
    echo "#SBATCH -N 1" >> "$filename"
    echo "#SBATCH --mem=64G" >> "$filename"
    echo "#SBATCH --ntasks-per-node=15" >> "$filename"
    echo "#SBATCH -o %x.o%j" >> "$filename"
    echo "#SBATCH -e %x.e%j" >> "$filename"
    echo "module load samtools" >> "$filename"
    
    sort -V $item.minibamlist.txt > $item.minibamlist.ordered.txt
    cat $item.minibamlist.ordered.txt | xargs -n 10 > $item.minibamlist.ordered.grouped.txt
    line_count=$(wc -l < $item.minibamlist.ordered.grouped.txt)
    echo "Number of lines in $item.minibamlist.ordered.grouped.txt: ${line_count}"
    seq -f "samtools merge -@ 8 ./merge.R1/merged_R1_%04g.bam" 1 ${line_count} > "$item.headers.txt"
    paste -d ' ' $item.headers.txt $item.minibamlist.ordered.grouped.txt >> "$filename"

    echo "R1 merge script created for $item"
 
    # Cleanup
    # rm "${item}.minibamlist.txt"  # Uncomment if you want to remove the file
    rm "${item}.minibamlist.ordered.txt"
    rm "${item}.minibamlist.ordered.grouped.txt"
    rm "$item.headers.txt"

    # Move script into directory
    mv "$filename" ${sample_directory}

    # Move to the merge directory for the current item
    cd ${sample_directory}
         
    filename2="${item}.mergeR2scripts.sh"
    echo "#! /bin/bash" > "$filename2"
    echo "# Bercin Cenik" >> "$filename2"
    echo "#SBATCH -A b1042" >> "$filename2"
    echo "#SBATCH -p genomics" >> "$filename2"
    echo "#SBATCH -t 24:00:00" >> "$filename2"
    echo "#SBATCH -N 1" >> "$filename2"
    echo "#SBATCH --mem=64G" >> "$filename2"
    echo "#SBATCH --ntasks-per-node=15" >> "$filename2"
    echo "#SBATCH -o %x.o%j" >> "$filename2"
    echo "#SBATCH -e %x.e%j" >> "$filename2"
    echo """ >> "$filename2"
    echo "module load samtools" >> "$filename2"
    echo """ >> "$filename2"

    cat $item.scripts.sh | tail -n +14 | cut -f 5 -d " " > $item.R1list.txt
    sort -V $item.R1list.txt > $item.R1list.ordered.txt
    cat $item.R1list.ordered.txt | xargs -n 10 > $item.R1list.ordered.grouped.txt
    line_count=$(wc -l < $item.R1list.ordered.grouped.txt) 
    echo "Number of lines in $item.R1list.ordered.grouped.txt: ${line_count}"
    seq -f "samtools merge -@ 8 ./merge.R2/merged_R2_%03g.bam" 1 ${line_count} > "$item.mergeR2.headers.txt"
    paste -d ' ' $item.mergeR2.headers.txt $item.R1list.ordered.grouped.txt >> "$filename2"
    
    echo "R2 merge script created for ${item}"

    # make the grouped R2 lists
    cat "${item}.mergeR2scripts.sh" | tail -n +14 | cut -f 5 -d " " > ${item}.R2list.txt
    sort -V ${item}.R2list.txt > ${item}.R2list.ordered.txt
    cat "${item}.R2list.ordered.txt" | xargs -n 10 > ${item}.R2list.ordered.grouped.txt

    # make the R3 merge scripts
    filename3="${item}.mergeR3scripts.sh"
    echo "#! /bin/bash" > "$filename3"
    echo "# Bercin Cenik" >> "$filename3"
    echo "#SBATCH -A b1042" >> "$filename3"
    echo "#SBATCH -p genomics" >> "$filename3"
    echo "#SBATCH -t 24:00:00" >> "$filename3"
    echo "#SBATCH -N 1" >> "$filename3"
    echo "#SBATCH --mem=64G" >> "$filename3"
    echo "#SBATCH --ntasks-per-node=15" >> "$filename3"
    echo "#SBATCH -o %x.o%j" >> "$filename3"
    echo "#SBATCH -e %x.e%j" >> "$filename3"
    echo "" >> "$filename3"
    echo "module load samtools" >> "$filename3"
    echo "" >> "$filename3"

    # Count lines from R2 grouped file and generate merge commands for R3
    line_count=$(wc -l < $item.R2list.ordered.grouped.txt)
    echo "Number of lines is ${line_count}"
    seq -f "samtools merge -@ 8 ./merge.R3/merged_R3_%03g.bam" 1 ${line_count} > "$item.mergeR3.headers.txt"
    paste -d ' ' $item.mergeR3.headers.txt $item.R2list.ordered.grouped.txt >> "$filename3"

    echo "R3 merge script created for ${item}"

    # make the grouped R3 lists
    cat "${item}.mergeR3scripts.sh" | tail -n +14 | cut -f 5 -d " " > ${item}.R3list.txt
    sort -V ${item}.R3list.txt > ${item}.R3list.ordered.txt
    cat "${item}.R3list.ordered.txt" | xargs -n 10 > ${item}.R3list.ordered.grouped.txt

    filename4="${item}.mergeR4scripts.sh"
    echo "#! /bin/bash" > "$filename4"
    echo "# Bercin Cenik" >> "$filename4"
    echo "#SBATCH -A b1042" >> "$filename4"
    echo "#SBATCH -p genomics" >> "$filename4"
    echo "#SBATCH -t 24:00:00" >> "$filename4"
    echo "#SBATCH -N 1" >> "$filename4"
    echo "#SBATCH --mem=64G" >> "$filename4"
    echo "#SBATCH --ntasks-per-node=15" >> "$filename4"
    echo "#SBATCH -o %x.o%j" >> "$filename4"
    echo "#SBATCH -e %x.e%j" >> "$filename4"
    echo "module load samtools" >> "$filename4"
    echo "" >> "$filename4"

    # Count lines from R3 grouped file and generate merge commands for R4
    line_count=$(wc -l < $item.R3list.ordered.grouped.txt)
    echo "Number of lines is ${line_count}"
    seq -f "samtools merge -@ 8 ./merge.R4/merged_R4_%03g.bam" 1 ${line_count} > "$item.mergeR4.headers.txt"
    paste -d ' ' $item.mergeR4.headers.txt $item.R3list.ordered.grouped.txt >> "$filename4"

    # Post-process and move the final bam into the working directory
    echo "" >> "$filename4"
    echo "final_bam_path=\"${sample_directory}/merge.R4/merged_R4_001.bam\"" >> "$filename4"
    echo "if [ -f \"\$final_bam_path\" ]; then" >> "$filename4"
    echo "    final_destination=\"${working_directory}/\$item.finalmerge.bam\"" >> "$filename4"
    echo "    # Check if the final file exists, if so, generate a new name" >> "$filename4"
    echo "    counter=1" >> "$filename4"
    echo "    while [ -f \"\$final_destination\" ]; do" >> "$filename4"
    echo "        final_destination=\"${working_directory}/\$item.finalmerge_\$counter.bam\"" >> "$filename4"
    echo "        counter=\$((counter+1))" >> "$filename4"
    echo "    done" >> "$filename4"
    echo "    mv \"\$final_bam_path\" \"\$final_destination\"" >> "$filename4"
    echo "    echo \"Final BAM file moved to \$final_destination\"" >> "$filename4"
    echo "else" >> "$filename4"
    echo "    echo \"File not found: \$final_bam_path\"" >> "$filename4"
    echo "fi" >> "$filename4"

    echo "R4 merge script with post-processing created for ${item}"

    echo "cleaning up.."
    rm *.txt

    ## submit to scheduler based on user input
    # Check the run_flag
    if [ "$run_flag" -eq 1 ]; then
        job1=$(sbatch "$item.scripts.sh" | awk '{print $NF}')
        echo "Job ID for ${item}.scripts.sh: $job1"

        job2=$(sbatch --dependency=afterok:$job1 "$item.mergeR2scripts.sh" | awk '{print $NF}')
        echo "Job ID for ${item}.mergeR2scripts.sh: $job2"    
        
        job3=$(sbatch --dependency=afterok:$job2 "$item.mergeR3scripts.sh" | awk '{print $NF}')
        echo "Job ID for ${item}.mergeR3scripts.sh: $job3"

        job4=$(sbatch --dependency=afterok:$job3 "$item.mergeR4scripts.sh" | awk '{print $NF}')
        echo "Job ID for ${item}.mergeR4scripts.sh: $job4"
    elif [ "$run_flag" -eq 0 ]; then
        echo "Scripts created for $item"
    else
        echo "Invalid -run value. Only 0 or 1 allowed."
        exit 1
    fi
    
    # Return to the parent directory
    cd ${working_directory}

    echo "Done with ${item}"
done < "$bampaths"

echo "finished running BA-BAM"
