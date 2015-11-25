##Primer Tool for Assembly using Yeast Protocol 

This is a simple program that reads in a file with multiple sequences in fasta
format and then finds the primers that are required to join these sequences
together using a yeast assembly protocol.

usage: ```find_primers.py [-h] [-oL OVERLAP_LENGTH] [-eL ELONGATION_LENGTH]  inputFile outputFile```


#### Positional arguments:
    inputFile
    outputFile

#### Optional arguments:
    -h, --help                show the help message and exit
    -oL, --overlap_length     OVERLAP_LENGTH    - this is the number of b.p of the primer that overlaps the adjacent segment        
    -eL, --elongation_length  ELONGATION_LENGTH - this is the nmber of b.p in the elongation part of the primer


###Dependency
primer3 (https://github.com/libnano/primer3-py)

You can use the get_primer3.sh script to download and install the primer3 module. This will be done inplace where this repo was cloned.
