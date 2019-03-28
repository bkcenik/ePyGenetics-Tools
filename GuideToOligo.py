
## GuideToOligo ##
## biopython wrapper script for converting your guide RNA sequences to order-ready oligos
## coded by Bercin Cenik

# create virtual environment, install BioPython on it
# here are the shell commands for that
    ## conda create --name BioPython
        # create virtual environment, can clone from existing python env
        # Note: I use 3.6.7
    ## conda activate BioPython
    # install biopython and/or update it
        ## pip install biopython
        ## pip install biopython --upgrade
    # you are now ready to use biopython

# to run this script:
# python GuideToOligo.py


## GuideToOligo Python Script ##

# import dependencies
import pandas as pd
from Bio.Seq import Seq

# make a list with your guides in the following format:
# Header: Name, Sequence
# input guide names and sequences
# do not include PAM in the guide
# save your file as a csv

# read in your guide sequences
prompt = input(f"Please input the path to your csv file: ")
guides = pd.read_csv(prompt, encoding="utf-8")
guides

# initialize empty list
column1 = []
column2 = []

# loop through input to create forward and reverse oligo sequences:
for x in guides["Sequence"]:
    Forward_oligo = "CACCG" + x
    column1.append(Forward_oligo)
    Reverse_oligo = "AAAC" + str(Seq(x).reverse_complement()) + "C"
    column2.append(Reverse_oligo)

# create new dataframe
guides2 = pd.DataFrame({
    "Name": guides["Name"],
    "Sequence": guides["Sequence"],
    "Forward Oligo": column1,
    "Reverse Oligo": column2
})

# print dataframe to console
guides2

# output to csv file
guides2.to_csv("output.csv")

