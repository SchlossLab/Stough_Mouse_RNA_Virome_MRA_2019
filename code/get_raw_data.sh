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

	echo "Downloading $accession ..."
	echo
	fasterq-dump --split-3 -O data/raw/raw/ "$accession"
	mv data/raw/raw/"$accession"_1.fastq data/raw/raw/"$accession"_forward.fastq
	mv data/raw/raw/"$accession"_2.fastq data/raw/raw/"$accession"_reverse.fastq

done

rename SRR5124216_ germ_free_ data/raw/raw/*.fastq
rename SRR5124214_ streptomycin_630_ data/raw/raw/*.fastq
rename SRR5124213_ clindamycin_630_ data/raw/raw/*.fastq
rename SRR5124212_ cefoperazone_630_ data/raw/raw/*.fastq
rename SRR6216945_ streptomycin_mock_ data/raw/raw/*.fastq
rename SRR6216939_ clindamycin_mock_ data/raw/raw/*.fastq
rename SRR6216941_ cefoperazone_mock_ data/raw/raw/*.fastq
