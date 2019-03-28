
# create virtual environment
    ## conda create --name BioPython
        # create virtual environment, can clone from existing python env
        # Note: I use 3.6.7
    ## conda activate BioPython
    # install biopython and/or update it
        ## pip install biopython
        ## pip install biopython --upgrade
    # you are now ready to use biopython

print("Hello, welcome to revcomp!")
print("-------------------------------------------------------------------")

from Bio.Seq import Seq
my_seq = Seq(input("please input your sequence "))
revcomp = my_seq.reverse_complement()

print("-------------------------------------------------------------------")
print("The reverse complement of your sequence is:")
print(revcomp)

print("-------------------------------------------------------------------")
print("Thank you for using revcomp ^_^")
