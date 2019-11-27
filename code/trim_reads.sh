#!/usr/bin/bash

####################################################################
# This script prepares raw short reads for assembly by removing    #
# read pairs that map to a collection of rRNA databases, removing, #
# lowquality reads, trimming sequencing adapters and low quality   #
# bases.                                                           #
# To run this script you must have Trimmomatic, BBMap, and         #   
# SortMeRNA installed.                                             #
####################################################################

### Set up environment

raw_path="data/raw/raw"
trimmed_path="data/raw/trimmed"
ref_fa_path="data/references/sortmerna_fa"
ref_db_path="data/references/sortmerna_db"

### This chunk checks for whether the trimmed data directory exists
# and creates it if it doesn't. This directory is necessary for the
# trimming step in the code below.

if [ -d $trimmed_path ]
then
	echo "Trimmed read folder already exists, continuing..."
	echo
else
	echo "Trimmed read folder doesn't exist, creating and continuing..."
	echo
	mkdir $trimmed_path
fi

### This chunk checks for the existence of the sortmeRNA reference
# databases, and creates them if they don't. 

if [ -d $ref_fa_path ]
then
    	echo "rRNA ref folder already exists, continuing..."
        echo
else
    	echo "rRNA ref folder doesn't exist, creating and continuing..."
        echo
	mkdir $ref_fa_path
fi

if [ -d $ref_db_path ]
then
    	echo "rRNA db folder already exists, continuing..."
        echo
else
    	echo "rRNA db folder doesn't exist, creating and continuing..."
        echo
	mkdir $ref_db_path
fi

### This chunk downloads the sortmeRNA rRNA databases into the directory
# created above.

wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/rfam-5.8s-database-id98.fasta -O $ref_fa_path/rfam-5.8s-database-id98.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/rfam-5s-database-id98.fasta -O $ref_fa_path/rfam-5s-database-id98.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-arc-16s-id95.fasta -O $ref_fa_path/silva-arc-16s-id95.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-arc-23s-id98.fasta -O $ref_fa_path/silva-arc-23s-id98.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-bac-16s-id90.fasta -O $ref_fa_path/silva-bac-16s-id90.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-bac-23s-id98.fasta -O $ref_fa_path/silva-bac-23s-id98.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-euk-18s-id95.fasta -O $ref_fa_path/silva-euk-18s-id95.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-euk-28s-id98.fasta -O $ref_fa_path/silva-euk-28s-id98.fasta

### This chunk creates indexed databases from the rRNA fasta files downloaded
# above and deposits them into the directory created above.

indexdb_rna --ref $ref_fa_path/silva-bac-16s-id90.fasta,$ref_db_path/silva-bac-16s-db:$ref_fa_path/silva-bac-23s-id98.fasta,$ref_db_path/silva-bac-23s-db:$ref_fa_path/silva-arc-16s-id95.fasta,$ref_db_path/silva-arc-16s-db:$ref_fa_path/silva-arc-23s-id98.fasta,$ref_db_path/silva-arc-23s-db:$ref_fa_path/silva-euk-18s-id95.fasta,$ref_db_path/silva-euk-18s-db:$ref_fa_path/silva-euk-28s-id98.fasta,$ref_db_path/silva-euk-28s:$ref_fa_path/rfam-5s-database-id98.fasta,$ref_db_path/rfam-5s-db:$ref_fa_path/rfam-5.8s-database-id98.fasta,$ref_db_path/rfam-5.8s-db

### This chunk loops through a list of sample names contained in
# data/process/samples.tsv and uses BBMap to join paired read files
# into a single interleaved fastq format, then uses sortmeRNA to
# remove read pairs that map to the rRNA databases created above.
# I then uses Trimmomatic to trim low quality bases and adapter 
# sequences from raw short read files. The loop then concatenates 
# unpaired reads produced by trimming into a single unpaired read file.

for sample in $(awk '{ print $2 }' data/process/samples.tsv); do
	
	reformat.sh in1=$raw_path/"$sample"_forward.fastq in2=$raw_path/"$sample"_reverse.fastq out=$raw_path/"$sample".fastq overwrite=true
	sortmerna --ref data/references/sortmerna_fa/silva-bac-16s-id90.fasta,$ref_db_path/silva-bac-16s-db:$ref_fa_path/silva-bac-23s-id98.fasta,$ref_db_path/silva-bac-23s-db:$ref_fa_path/silva-arc-16s-id95.fasta,$ref_db_path/silva-arc-16s-db:$ref_fa_path/silva-arc-23s-id98.fasta,$ref_db_path/silva-arc-23s-db:$ref_fa_path/silva-euk-18s-id95.fasta,$ref_db_path/silva-euk-18s-db:$ref_fa_path/silva-euk-28s-id98.fasta,$ref_db_path/silva-euk-28s:$ref_fa_path/rfam-5s-database-id98.fasta,$ref_db_path/rfam-5s-db:$ref_fa_path/rfam-5.8s-database-id98.fasta,$ref_db_path/rfam-5.8s-db --reads $raw_path/"$sample".fastq --aligned $raw_path/"$sample"_rrna --other $raw_path/"$sample"_clean --paired_in --fastx -a 12
	reformat.sh in=$raw_path/"$sample"_clean.fastq out1=$raw_path/"$sample"_forward.fastq out2=$raw_path/"$sample"_reverse.fastq overwrite=true
	rm $raw_path/"$sample"_clean.fastq $raw_path/"$sample"_rrna.fastq
        trimmomatic PE -threads 12 $raw_path/"$sample"_forward.fastq $raw_path/"$sample"_reverse.fastq $trimmed_path/"$sample"_forward_paired.fastq $trimmed_path/"$sample"_forward_unpaired.fastq $trimmed_path/"$sample"_reverse_paired.fastq $trimmed_path/"$sample"_reverse_unpaired.fastq LEADING:10 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:50 ILLUMINACLIP:data/references/TruSeq3-PE.fa:2:30:10
	cat $trimmed_path/"$sample"_forward_unpaired.fastq $trimmed_path/"$sample"_reverse_unpaired.fastq > $trimmed_path/"$sample"_unpaired.fastq
	rm $trimmed_path/"$sample"_forward_unpaired.fastq $trimmed_path/"$sample"_reverse_unpaired.fastq

done
