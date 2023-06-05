#!/bin/bash
# make minibam file compression scripts: maketar.sh
# written by: Bercin Cenik

#### instructions
# create a CSV file with CHACHA number, e.g., CHACHAx.csv, save it in /data/data_offload/CHACHAx/
# header of CSV file should be: pool,firstbarcode,lastbarcode
# then enter the pool number and barcodes in corresponding rows. eg. 1,01,04 etc.
# in CHACHAx, run the script eg: sh maketar.sh
# it will ask you for the CHACHA number eg. x, and then the path to the csv file (type CHACHAx.csv)
# a tar merging script will be created eg. pool1.tar.sh
# each pool will have its own script.
# run these scripts individually, for example: sh pool1.tar.sh
# this will compress the minibam files as such: pool<x>.barcode<y>.bam.tar.gz

# Prompt the user for the path to the CSV file
read -p "Enter the path to the CSV file: " csv_file

# Check if the CSV file exists
if [ ! -f "$csv_file" ]; then
    echo "Error: CSV file not found."
    exit 1
fi

# Prompt the user for the CHACHA number
read -p "Enter the CHACHA number: " chacha_number

# Read the CSV file and process each line
tail -n +2 "$csv_file" | while IFS=',' read -r pool first last; do
    # Generate the pool tar file path
    pool_tar_file="pool${pool}.tar.sh"

    # Create the pool tar script
    echo "#!/bin/bash" > "$pool_tar_file"
    echo " " >> "$pool_tar_file"
    echo "# Pool ${pool} tar compression script" >> "$pool_tar_file"
    echo " " >> "$pool_tar_file"

    for barcode in $(seq -w "$first" "$last"); do
        # Generate the individual tar file path
        tar_file="pool${pool}.barcode${barcode}.bam.tar.gz"

        # Append the tar compression command to the pool tar script
        echo "tar cvzf ${tar_file} /data/data_offload/CHACHA${chacha_number}/CHACHA${chacha_number}_P${pool}/*/bam_pass/barcode${barcode}/*.bam" >> "$pool_tar_file"

        # Delete the individual barcode tar script
       # rm "pool${pool}.barcode${barcode}.bam.tar.gz.sh"
    done

    echo "echo \"Finished compressing all barcodes for Pool ${pool}\"" >> "$pool_tar_file"
done

cat pool*.sh > CHACHA5.tarminibams.sh
rm pool*.sh