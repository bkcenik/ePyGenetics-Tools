#!/bin/bash
# make merged fastq scripts: make.merged.fastq.sh
# written by: Bercin Cenik

#### instructions
# create a CSV file with CHACHA number, e.g., CHACHAx.csv, save it in /data/data_offload/CHACHAx/
# header of CSV file should be: pool,firstbarcode,lastbarcode
# then enter the pool number and barcodes in corresponding rows. eg. 1,01,04 etc.
# in CHACHAx, run the script: sh make.merged.fastq.sh
# it will ask you for the CHACHA number, e.g., x, and then the path to the csv file (type CHACHAx.csv)
# a fastq merging script will be created, with each pool having its own script.

# Prompt for CHACHA number
read -p "Enter CHACHA number (x): " chacha

# Prompt for CSV file path
read -p "Enter the path of the CSV file: " csv_file

# Create the CHACHAx.fastq.merge.sh script
merge_script_filename="CHACHA${chacha}.fastq.merge.sh"
echo "#!/bin/bash" > "$merge_script_filename"

# Read the CSV file and process each line
while IFS=',' read -r pool start end
do
  # Skip processing the first line (header)
  if [ "$pool" = "pool" ]; then
    continue
  fi

  # Generate the desired command for each barcode in the range
  for barcode in $(seq -w "$start" "$end")
  do
    echo "for fastq in \$(ls /data/data_offload/CHACHA${chacha}/CHACHA_P${pool}/*/fastq_pass/barcode${barcode}/fastq_pass/*.fastq.gz); do" >> "$merge_script_filename"
    echo "  cat \"\$fastq\" >> pool${pool}.barcode${barcode}.merge.fastq.gz" >> "$merge_script_filename"
    echo "done" >> "$merge_script_filename"
  done

  echo "Script commands for pool${pool} have been added to '$merge_script_filename'."
done < "$csv_file"

# Make the merge script file executable
chmod +x "$merge_script_filename"

echo "The CHACHA${chacha}.fastq.merge.sh script has been created."
