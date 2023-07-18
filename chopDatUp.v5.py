## Chop Dat Up
## author: Bercin Cenik
## this script is for splitting fasta files >2kb into smaller chunks so they can be submitted to CRISPick

import os
import argparse
import logging
import sys

print("starting Chop-dat-Up")

def read_fasta(file_path):
    with open(file_path, "r") as file:
        header = file.readline().strip()
        sequence = "".join(line.strip() for line in file.readlines())
    return header, sequence

def write_fasta(file_path, header, sequence, batch_num):
    with open(file_path, "w") as file:
        new_header = f"{header}.chop{batch_num}"
        file.write(f"{new_header}\n")
        n = 60  # Number of characters per line in the sequence
        for i in range(0, len(sequence), n):
            file.write(sequence[i:i+n] + "\n")

def chop_sequence(sequence, batch_size):
    for i in range(0, len(sequence), batch_size):
        yield sequence[i:i + batch_size]

def chop_fasta(input_file, output_dir, max_batch_size=2000):
    header, sequence = read_fasta(input_file)
    seq_length = len(sequence)

    if seq_length <= max_batch_size:
        output_file = os.path.join(output_dir, os.path.basename(input_file))
        write_fasta(output_file, header, sequence, batch_num="")
        logging.info(f"{output_file} has {seq_length} bases. Not chopping.")
    else:
        logging.info(f"{input_file} has {seq_length} bases. Chopping.")

        batch_num = 1
        total_files = 0
        for batch_seq in chop_sequence(sequence, max_batch_size):
            output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}.chop{batch_num}.fa")
            write_fasta(output_file, header, batch_seq, batch_num=batch_num)
            logging.info(f"Chopped into {output_file}.")
            batch_num += 1
            total_files += 1

        logging.info(f"{input_file} was split into {total_files} files.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Chop up long fasta files into smaller batches.")
    parser.add_argument("input_dir", help="Path to the input directory containing the fasta files.")
    parser.add_argument("output_dir", help="Path to the output directory where the chopped fasta files will be saved.")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir

    if not os.path.exists(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist.")
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    logging.basicConfig(filename='chop_fasta.log', level=logging.INFO,
                        format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    for filename in os.listdir(input_dir):
        if filename.endswith(".fa"):
            input_file = os.path.join(input_dir, filename)
            chop_fasta(input_file, output_dir)

    # Move the chopped files to subdirectories based on the number of files each entry was split into
    for file in os.listdir(output_dir):
        if file.endswith(".fa"):
            num_splits = int(file.split(".chop")[1][0])  # Extract the number of splits from the filename
            subfolder_path = os.path.join(output_dir, f"split.{num_splits}")
            if not os.path.exists(subfolder_path):
                os.makedirs(subfolder_path)
            os.rename(os.path.join(output_dir, file), os.path.join(subfolder_path, file))

print("Chop-dat-Up complete. Please check your outputs.")