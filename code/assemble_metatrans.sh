#!/usr/bin/bash

####################################################################
### This script assembles RNAseq reads into contigs and scaffolds  #
# To run this script you must have SPAdes installed.               #
####################################################################


### This chunk checks for whether the assembly data directory exists
# and creates it if it doesn't. This directory is necessary for the
# trimming step in the code below.

if [ -d "data/raw/assembly/" ]
then
    	echo "Assembly folder already exists, continuing..."
        echo
else
    	echo "Assembly folder doesn't exist, creating and continuing..."
        echo
	mkdir data/raw/assembly
fi


### This chunk creates a temporary directory for use by SPAdes during assembly 
# and saves it in the $TEMP variable. The if statement exits the script with 
# an error if the temp directory has not been successfully created.

TEMP=$(mktemp -d)
if [ $? -ne 0 ]
then
    	echo "can't create $TEMP ... exit"
        exit 1
fi


### This chunk iterates over a list of treatment groups contained in the
# data/process/samples.tsv file and assembles their respective short
# read sequence files using SPAdes in RNA mode with a maximum RAM limit
# of 900GB with 32 threads which are then stored in separate directories
# within data/raw/assembly/. The last line deletes the $TEMP directory
# created above.

for treatment in $(awk '{ print $2 }' data/process/samples.tsv); do

	mkdir data/raw/assembly/"$treatment"
	spades.py --rna -o data/raw/assembly/"$treatment"/ -1 data/raw/trimmed/"$treatment"_forward_paired.fastq -2 data/raw/trimmed/"$treatment"_reverse_paired.fastq -s data/raw/trimmed/"$treatment"_unpaired.fastq -t 32 -m 900 --tmp-dir $TEMP

done

rm -rf $TEMP
