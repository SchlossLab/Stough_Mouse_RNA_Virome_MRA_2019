#!/usr/bin/bash

####################################################################
# This script prepares raw short reads for assembly and mapping.   #
# To run this script you must have Trimmomatic installed.          #
####################################################################


### This chunk checks for whether the trimmed data directory exists
# and creates it if it doesn't. This directory is necessary for the
# trimming step in the code below.

if [ -d "data/raw/trimmed/" ]
then
	echo "Trimmed read folder already exists, continuing..."
	echo
else
	echo "Trimmed read folder doesn't exist, creating and continuing..."
	echo
	mkdir data/raw/trimmed
fi

### This chunk loops through a list of sample names contained in
# data/process/samples.tsv and uses Trimmomatic to trim low quality 
# bases and adapter sequences from raw short read files. The loop
# then concatenates unpaired reads produced by trimming into a
# single unpaired read file.

for sample in $(awk '{ print $2 }' data/process/samples.tsv); do

        trimmomatic PE -threads 12 data/raw/raw/"$sample"_forward.fastq data/raw/raw/"$sample"_reverse.fastq data/raw/trimmed/"$sample"_forward_paired.fastq data/raw/trimmed/"$sample"_forward_unpaired.fastq data/raw/trimmed/"$sample"_reverse_paired.fastq data/raw/trimmed/"$sample"_reverse_unpaired.fastq LEADING:10 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:50 ILLUMINACLIP:data/references/TruSeq3-PE.fa:2:30:10
	cat data/raw/trimmed/"$sample"_forward_unpaired.fastq data/raw/trimmed/"$sample"_reverse_unpaired.fastq > data/raw/trimmed/"$sample"_unpaired.fastq
	rm data/raw/trimmed/"$sample"_forward_unpaired.fastq data/raw/trimmed/"$sample"_reverse_unpaired.fastq

done


