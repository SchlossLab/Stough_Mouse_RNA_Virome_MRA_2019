#!/usr/bin/bash

###########################################################################
# This script processes maps trimmed short reads to scaffolds for         # 
# quantification of relative abundance.                                   #
# A working installation of Bowtie2 is required to                        #
# to run this script.                                                     #
###########################################################################

### This chunk uses the build function from Bowtie 2 to make index files
# from the whole scaffold file for mapping short reads to scaffolds.

bowtie2-build data/raw/scaffolds/all_scaffolds.fasta data/mapping/index/all_scaffolds 

### This chunck runs read mappings, converts SAM files to BAM format, 
# and sorts the BAM files.

for $treatment in $(awk '{ print $2 }' data/process/samples.txt); do

	bowtie2 -x data/raw/mapping/index/all_scaffolds -1 data/raw/trimmed/"$treatment"_forward_paired.fastq\
	-2 data/raw/trimmed/"$treatment"_reverse_paired.fastq -U data/raw/trimmed/"$treatment"_unpaired.fastq\
	-S data/raw/mapping/"$treatmenr".sam --end-to-end --sensitive -p 12 -t
	samtools view -S -b data/raw/mapping/"$treatment".sam > data/raw/mapping/"$treatment".bam
	samtools sort data/raw/mapping/"$treatment".bam -o data/raw/mapping/"$treatment"_sorted.bam
	rm data/raw/mapping/"$treatment".sam

done

