# making merged fastq files and compressed minibam files, for NGS techs

Created by: Bercin Cenik
Date: June 5, 2023 2:21 PM

Overview: this is a set of instructions for creating merged fastq files and compressing minibam files for ONT long-read sequencing runs, written for NGS technicians in the Shilatifard lab.

Ensure that when creating the experiments, you use the following directory structure:

CHACHAx

CHACHA_Pm

CHACHA_Pn

…

where x is the CHACHA# and m and n are the barcode number. Always save experiments in data_offload.

Example: as part of CHACHA5, there are 5 pools, each containing 4 barcoded samples, being run on 5 separate flow cells. This folder structure would be:

- CHACHA5
    - CHACHA5_P1
        - fastq_pass
            - barcode01
                
                `*.fastq.gz` 
                
            - barcode02
                
                `*.fastq.gz` 
                
            - barcode03
                
                `*.fastq.gz` 
                
            - barcode04
                
                `*.fastq.gz` 
                
        - bam_pass
            - barcode01
                
                `*.bam` 
                
            - barcode02
                
                `*.bam` 
                
            - barcode03
                
                `*.bam` 
                
            - barcode04
                
                `*.bam` 
                
    - CHACHA5_P2
        
        fastq_pass..
        
        bam_pass..
        

Create a .csv file in excel or a text editor, in the following format:

```
pool,firstbarcode,lastbarcode
1,01,04
2,05,08
3,09,12
4,13,16
5,17,20
```

This looks like a table similar to the following:

| pool | firstbarcode | lastbarcode |
| --- | --- | --- |
| 1 | 01 | 04 |
| 2 | 05 | 08 |
| 3 | 09 | 12 |
| 4 | 13 | 16 |
| 5 | 17 | 20 |

Save this as CHACHAx.csv and move it to the experiment directory.

Copy `make.merged.fastq.sh` and `make.tar.sh` into the experiment directory as well.

### Merge and transfer merged fastq files to b1042

`make.merged.fastq.sh`

```bash
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
```

Run this script:

```bash
sh make.merged.fastq.sh
```

It will prompt the user for the CHACHA# and the path to the CSV file with pool and barcode numbers, and generate the merge script:

```bash
>Enter CHACHA number (x): 5
>Enter the path of the CSV file: CHACHA5.csv
>Script commands for pool1 have been added to 'CHACHA5.fastq.merge.sh'.
>Script commands for pool2 have been added to 'CHACHA5.fastq.merge.sh'.
>Script commands for pool3 have been added to 'CHACHA5.fastq.merge.sh'.
>Script commands for pool4 have been added to 'CHACHA5.fastq.merge.sh'.
>Script commands for pool5 have been added to 'CHACHA5.fastq.merge.sh'.
>The CHACHA5.fastq.merge.sh script has been created.
```

These are the contents of the `CHACHAx.fastq.merge.sh` script:

```bash
#!/bin/bash
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P1/*/fastq_pass/barcode01/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool1.barcode01.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P1/*/fastq_pass/barcode02/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool1.barcode02.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P1/*/fastq_pass/barcode03/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool1.barcode03.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P1/*/fastq_pass/barcode04/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool1.barcode04.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P2/*/fastq_pass/barcode05/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool2.barcode05.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P2/*/fastq_pass/barcode06/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool2.barcode06.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P2/*/fastq_pass/barcode07/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool2.barcode07.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P2/*/fastq_pass/barcode08/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool2.barcode08.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P3/*/fastq_pass/barcode09/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool3.barcode09.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P3/*/fastq_pass/barcode10/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool3.barcode10.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P3/*/fastq_pass/barcode11/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool3.barcode11.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P3/*/fastq_pass/barcode12/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool3.barcode12.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P4/*/fastq_pass/barcode13/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool4.barcode13.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P4/*/fastq_pass/barcode14/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool4.barcode14.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P4/*/fastq_pass/barcode15/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool4.barcode15.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P4/*/fastq_pass/barcode16/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool4.barcode16.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P5/*/fastq_pass/barcode17/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool5.barcode17.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P5/*/fastq_pass/barcode18/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool5.barcode18.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P5/*/fastq_pass/barcode19/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool5.barcode19.merge.fastq.gz
done
for fastq in $(ls /data/data_offload/CHACHA5/CHACHA_P5/*/fastq_pass/barcode20/fastq_pass/*.fastq.gz); do
  cat "$fastq" >> pool5.barcode20.merge.fastq.gz
done
```

Running this script (`sh CHACHAx.fastq.merge.sh`) will output the merged fastq files.

While still in the experiment directory, make a merged fastq directory and move all the files there.

```bash
cd /data/data_offload/CHACHAx
mkdir fastq.merge
cd fastq.merge
mv ../*fastq.gz . 
```

Following this, login to Quest and create the directories for the experiments in b1042, then log out of Quest.

```bash
ssh -X <netID>@quest.it.northwestern.edu
# enter password when prompted
```

```bash
mkdir /projects/b1042/Shilatifard/PromethION/CHACHAx
cd CHACHAx
mkdir ./bam ./minibam ./tars ./scripts ./QC.reports ./seq.summary ./fastq

exit
```

Transfer the files to b1042. This will take a while:

```bash
cd /data/data_offload/CHACHAx/
scp -r ./fastq.merge <netID>@quest.it.northwestern.edu:/projects/b1042/Shilatifard/PromethION/CHACHAx/fastq
# enter password when prompted.
```

### Compress minibam files and transfer to b1042:

`make.tar.sh`

```bash
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
```

Run `make.tar.sh`:

```bash
sh make.tar.sh
```

This will create `CHACHA5.tarminibams.sh`:

```bash
!/bin/bash
 
# Pool 1 tar compression script
 
tar cvzf pool1.barcode01.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P1/*/bam_pass/barcode01/*.bam
tar cvzf pool1.barcode02.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P1/*/bam_pass/barcode02/*.bam
tar cvzf pool1.barcode03.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P1/*/bam_pass/barcode03/*.bam
tar cvzf pool1.barcode04.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P1/*/bam_pass/barcode04/*.bam
echo "Finished compressing all barcodes for Pool 1"
#!/bin/bash
 
# Pool 2 tar compression script
 
tar cvzf pool2.barcode05.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P2/*/bam_pass/barcode05/*.bam
tar cvzf pool2.barcode06.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P2/*/bam_pass/barcode06/*.bam
tar cvzf pool2.barcode07.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P2/*/bam_pass/barcode07/*.bam
tar cvzf pool2.barcode08.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P2/*/bam_pass/barcode08/*.bam
echo "Finished compressing all barcodes for Pool 2"
#!/bin/bash
 
# Pool 3 tar compression script
 
tar cvzf pool3.barcode09.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P3/*/bam_pass/barcode09/*.bam
tar cvzf pool3.barcode10.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P3/*/bam_pass/barcode10/*.bam
tar cvzf pool3.barcode11.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P3/*/bam_pass/barcode11/*.bam
tar cvzf pool3.barcode12.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P3/*/bam_pass/barcode12/*.bam
echo "Finished compressing all barcodes for Pool 3"
#!/bin/bash
 
# Pool 4 tar compression script
 
tar cvzf pool4.barcode13.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P4/*/bam_pass/barcode13/*.bam
tar cvzf pool4.barcode14.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P4/*/bam_pass/barcode14/*.bam
tar cvzf pool4.barcode15.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P4/*/bam_pass/barcode15/*.bam
tar cvzf pool4.barcode16.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P4/*/bam_pass/barcode16/*.bam
echo "Finished compressing all barcodes for Pool 4"
#!/bin/bash
 
# Pool 5 tar compression script
 
tar cvzf pool5.barcode17.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P5/*/bam_pass/barcode17/*.bam
tar cvzf pool5.barcode18.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P5/*/bam_pass/barcode18/*.bam
tar cvzf pool5.barcode19.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P5/*/bam_pass/barcode19/*.bam
tar cvzf pool5.barcode20.bam.tar.gz /data/data_offload/CHACHA5/CHACHA5_P5/*/bam_pass/barcode20/*.bam
echo "Finished compressing all barcodes for Pool 5"
```

Run `CHACHA5.tarminibams.sh`:

```bash
sh CHACHA5.tarminibams.sh
```

When finished, create and move compresed bams into the bam.tars directory.

```bash
cd /data/data_offload/CHACHAx
mkdir bam.tars
cd bam.tars
mv ../*bam.tar.gz *
```

Transfer into the minibam directory using scp, this will take a while:

```bash
cd /data/data_offload/CHACHAx/
scp -r ./bam.tars <netID>@quest.it.northwestern.edu:/projects/b1042/Shilatifard/PromethION/CHACHAx/minibam
# enter password when prompted.
```