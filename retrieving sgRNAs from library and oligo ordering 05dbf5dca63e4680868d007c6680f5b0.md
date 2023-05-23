# retrieving sgRNAs from library and oligo ordering

Experiment: CRISPR
Last edited: May 9, 2023 6:48 PM
Owner: Bercin Cenik
Type: csv, python, shell
Works?: No

### Notes and debugging

add notes about debugging here.

### Usage

1. Make an excel sheet with the following columns:

#1: geneID

#2: 20bp seed sequence

#3: ”CACCG” for seeds not beginning with a “G”, “CACC” for seeds beginning with a “G”

#4: merge of columns 3&2

#5: “AAAC”

#6: reverse complement of seed sequence

#7: ”C” for seeds not beginning with a “G”, blank for seeds beginning with a “G”

#8: merge of columns 5&6&7

Next, populate these columns.

1. Column#2: Get your sgRNAs from the Brunello library. Use my “grep sgRNA from list” protocol for this. For example:
    
    ```python
    cat broadgpp-brunello-library-corrected.txt | grep -iw 'ZC3H13\|KIAA1429\|METTL3\|METTL14\|ELP1\|FOS\|USP49\|SNRNP25\|OGT\|SRF\|AKIRIN2\|MAPK1\|DDX49\|RIOK2' > methylation.sgRNAs.txt
    ```
    
    This is assuming you are in a directory that contains the library (”broadgpp-brunello-library-corrected.txt”). The names of genes need to be separated by “\|”, this will pipe the output to the subsequent gene, first escaping the pipe character. 
    
    The text file will look like this:
    
    ```python
    AKIRIN2_41635_CTGGTGAAGCTGGTCCACTG	CTGGTGAAGCTGGTCCACTG	AKIRIN2
    AKIRIN2_41636_GCCGCAGAAGTATCTCCGAA	GCCGCAGAAGTATCTCCGAA	AKIRIN2
    AKIRIN2_41637_GGGGACGCCGGGCTCAACAG	GGGGACGCCGGGCTCAACAG	AKIRIN2
    AKIRIN2_41638_TGTGCATCAGAAGTACAACA	TGTGCATCAGAAGTACAACA	AKIRIN2
    
    #The third column of these include the seed sequence, 
    #which will be pasted into the new excel sheet’s second column (#2: 20bp seed sequence)
    ```
    
2. Column#3: This is the 5’ flank of the duplex oligo. If the sgRNA seed sequence begins with a G, type “CACC” here. If it doesn’t, you will need to add a G manually to create a GC-clamp to increase the efficiency of the oligo. For seeds not beginning with a G, type “CACCG” here.
3. Column#4: the 5’ oligo, merge the 3th and 2nd columns by using & on excel, in that order.
4. Column#5: 3’ flank. This is always a “CAAA”, but since we are writing 5’ to 3’, type “AAAC” in this column for every entry.
5. Column #6: ”revcomp of sgRNA”. These are the reverse complements of the seed sequence. Briefly:
    
    ```bash
    original.            ACGTACGT
    reverse.             TGCATGCA
    complement.          TGCATGCA
    reverse complement.  ACGTACGT
    ```
    
    Run the “## for getting reverse complements” script shown below. This assumes you have a list of sgRNA seed sequences ("sgRNA_list.csv”) in the directory you are working in, and outputs a reverse complement of those (”rc_sgRNA_list.csv”). Paste the contents of rc_sgRNA_list.csv into column#6 (”revcomp of sgRNA”).
    
6. Column #7: clamp. This is a “C” if the seed does not begin with a G, blank if it does (since we are manually adding the GC if it’s absent, etc).
7. Column #8: the 3’ oligo, merge the 5th, 6th and 7th columns with the & operator on excel.

1. The final sheet will look like this:
    
    ```python
    "gene","sgRNA sequence","5'flank","5'oligo for ordering","3'flank","revcomp of sgRNA","clamp","3'oligo for ordering"
    AKIRIN2_sg1,CTGGTGAAGCTGGTCCACTG,CACCG,CACCGCTGGTGAAGCTGGTCCACTG,AAAC,CAGTGGACCAGCTTCACCAG,C,AAACCAGTGGACCAGCTTCACCAGC
    AKIRIN2_sg2,GCCGCAGAAGTATCTCCGAA,CACC,CACCGCCGCAGAAGTATCTCCGAA,AAAC,TTCGGAGATACTTCTGCGGC,,AAACTTCGGAGATACTTCTGCGGC
    AKIRIN2_sg3,GGGGACGCCGGGCTCAACAG,CACC,CACCGGGGACGCCGGGCTCAACAG,AAAC,CTGTTGAGCCCGGCGTCCCC,,AAACCTGTTGAGCCCGGCGTCCCC
    AKIRIN2_sg4,TGTGCATCAGAAGTACAACA,CACCG,CACCGTGTGCATCAGAAGTACAACA,AAAC,TGTTGTACTTCTGATGCACA,C,AAACTGTTGTACTTCTGATGCACAC
    ```
    
    The 4th and 8th columns are the sequences of the oligos which will be ordered.
    

1. Run the script “##for reformatting the sgRNA oligos to an order-ready layout” below. This assumes there is an input.csv in your working directory that has the above 8-column format for all guides, and outputs the sequences in tandem (Forward followed by Reverse for each oligo). The list should look like this:
    
    ```python
    AKIRIN2_sg1_F	CACCGCTGGTGAAGCTGGTCCACTG
    AKIRIN2_sg1_R	AAACCAGTGGACCAGCTTCACCAGC
    AKIRIN2_sg2_F	CACCGCCGCAGAAGTATCTCCGAA
    AKIRIN2_sg2_R	AAACTTCGGAGATACTTCTGCGGC
    AKIRIN2_sg3_F	CACCGGGGACGCCGGGCTCAACAG
    AKIRIN2_sg3_R	AAACCTGTTGAGCCCGGCGTCCCC
    AKIRIN2_sg4_F	CACCGTGTGCATCAGAAGTACAACA
    AKIRIN2_sg4_R	AAACTGTTGTACTTCTGATGCACAC
    ```
    
2. Order oligos from IDT (25nm, STD).
3. Important note: make sure non-target controls are also in your list of oligos to order.

---

### Scripts

```python
## for getting reverse complements

import csv

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq = seq.lstrip('\ufeff')  # remove UTF-8 BOM character
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join([complement[base] for base in rev_seq])
    return rev_comp_seq

with open('sgRNA_list.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    sgRNA_list = [row[0] for row in reader]

rc_sgRNA_list = [reverse_complement(sgRNA) for sgRNA in sgRNA_list]

with open('rc_sgRNA_list.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for sgRNA in rc_sgRNA_list:
        writer.writerow([sgRNA])
```

```python
## for only getting reverse sequences:

import csv

def reverse_sequence(seq):
    return seq[::-1]

# Read sgRNA sequences from input CSV file
with open('three_five.csv', 'r') as f:
    reader = csv.reader(f)
    sgRNA_list = [row[0] for row in reader]

# Compute reverse of each sequence
reverse_sgRNA_list = [reverse_sequence(sgRNA) for sgRNA in sgRNA_list]

# Write reversed sequences to output CSV file
with open('reverse_sgRNAs.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    for sgRNA in reverse_sgRNA_list:
        writer.writerow([sgRNA])
```

```python
## for getting reverse complement of the nontarget controls:

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join([complement[base] for base in rev_seq])
    return rev_comp_seq
sgRNA_list = ['AAAAAGCTTCCGCCTGATGG',
              'AAAACAGGACGATGTGCGGC',
              'AAAACATCGACCGAAAGCGT',
              'AAAATAGCAGTAAACTCAAC']

rc_sgRNA_list = [reverse_complement(sgRNA) for sgRNA in sgRNA_list]

print(rc_sgRNA_list)
```

```python
##for reformatting the sgRNA oligos to an order-ready layout
import csv

with open('input.csv') as infile, open('output2.csv', 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile, delimiter='\t')
    for row in reader:
        for i in range(0, len(row), 2):
            writer.writerow([row[i], row[i+1]])
```

### Related files

[https://www.notion.so](https://www.notion.so)

[https://www.notion.so](https://www.notion.so)