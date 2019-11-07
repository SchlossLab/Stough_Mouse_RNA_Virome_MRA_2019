#!usr/bin/bash

###########################################################################
# This script processes assembled scaffolds, adding treatment groups to   #
# fasta sequence headers, removing short contigs, removing line breaks    #
# from fasta sequences, checks for circular contigs.                      #
# Working installations of BBMap and CContigs is required to              #
# to run this script.                                                     #
###########################################################################


### This chunk checks for whether the scaffold directory exists
# and creates it if it doesn't. This directory is necessary for the
# steps in the code below.

if [ -d "data/raw/scaffolds/" ]
then
    	echo "Scaffolds read folder already exists, continuing..."
        echo
else
    	echo "Scaffolds read folder doesn't exist, creating and continuing..."
        echo
        mkdir data/raw/scaffolds
fi

### This chunk adds filenames to fasta headers, concatenates scaffolds into a 
# single file, and removes scaffolds under 1000 bases long

echo "Removing contigs under 1kb, adding treatment groups to fasta headers, and concatenating contigs into single file"
echo

for treatment in $(awk '{ print $2 }' data/process/samples.tsv); do

        reformat.sh minlength=1000 overwrite=t in=data/raw/assembly/"$treatment"/transcripts.fasta out=data/raw/scaffolds/"$treatment"_scaffolds.fasta
        sed 's/^>/>'"$treatment"'_/g' data/raw/scaffolds/"$treatment"_scaffolds.fasta >> data/raw/scaffolds/all_long_scaffolds.fasta

done


### This chunk removes line breaks from the sequences for easy extraction
# later using grep

echo "Removing line breaks from the scaffold sequences and creating whole scaffold file"
echo
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' data/raw/scaffolds/all_long_scaffolds.fasta > data/raw/scaffolds/all_scaffolds.fasta


### This chunk checks for whether the cir_contigs directory exists
# and creates it if it doesn't. This directory is necessary for the
# steps in the code below.

if [ -d "data/process/cir_contigs" ]
then
    	echo "Circular contig folder already exists, continuing..."
        echo
else
    	echo "Circular folder doesn't exist, creating and continuing..."
        echo
	mkdir data/process/cir_contigs
fi

### This chunk calls on the ccontigs software to align the ends of
# scaffolds to each other to see whether they may be circular 
# molecules. The output is a tab-delimited table that contains 
# the scaffold ID.

julia $HOME/bin/ccontigs/ccontigs.jl -i data/raw/scaffolds/all_scaffolds.fasta -o data/process/cir_contigs/ccontigs_out.tsv

### This chunk extracts Circular contigs from whole scaffold file
# into a circular scaffold file using the output table from ccontigs 

for ccontig in $(awk '{ print $1 }' data/process/cir_contigs/ccontigs_out.tsv); do

	grep -A 1 "$ccontig" data/raw/scaffolds/all_scaffolds.fasta >> data/process/cir_contigs/cir_scaffolds.fasta

done

