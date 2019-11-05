#!/usr/bin/bash

#############################################################################
### This script downloads raw sequencing data into the data/raw/ directory. #
# To run this script you must have the SRA toolkit installed.               #
#############################################################################


### This chunk checks for whether the raw data directory exists
# and creates it if it doesn't.	This directory is necessary for	the
# trimming step	in the code below.

if [ -d	"data/raw/raw/" ]
then
        echo "Raw read folder already exists, continuing..."
        echo
else
        echo "Raw read folder doesn't exist, creating and continuing..."
        echo
        mkdir data/raw/raw 
fi


### Download Raw Sequences from SRA
# This chunk loops through a list of SRA accession numbers contained in 
# data/process/samples.tsv and uses fastq-dump from the SRA
# toolkit to download each of the SRA records as separate forward
# and reverse fastq files in the data/raw/raw/ directory. 

for accession in $(awk '{ print $1 }' data/process/samples.tsv); do

	for treatment in $(awk '{ print $2 }' data/process/samples.tsv); do

		echo "Downloading $accession ..."
		echo
		fasterq-dump --split-3 -o data/raw/raw/"$treatment" "$accession"
		echo
		echo "Renaming $accession.fastq to $treatment.fastq"
		echo
		mv data/raw/raw/"$treatment"_1.fastq data/raw/raw/"$treatment"_forward.fastq
		mv data/raw/raw/"$treatment"_2.fastq data/raw/raw/"$treatment"_reverse.fastq

	done

done
