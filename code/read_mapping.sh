#!/usr/bin/bash

###########################################################################
# This script processes assembled scaffolds, adding treatment groups to   #
# fasta sequence headers, removing short contigs, removing line breaks    #
# from fasta sequences, checks for circular contigs, and maps trimmed     #
# RNAseq short reads to scaffolds for quantification of relative          #
# abundance.                                                              #
# Working installations of Bowtie2, BBMap, and CContigs is required to    #
# to run this script.                                                     #
###########################################################################

### This chunk 

if [ -e $scaffold_path/all_raw_scaffolds.fasta ]
then
    	echo "Deleting old temp contig file"
        rm $scaffold_path/all_raw_scaffolds.fasta
else
    	echo "No old contig files to delete"
fi

echo

### Add filenames to fasta headers and concatenate contigs into single file

echo "Adding filenames to fasta headers and concatenating contigs into single file"
echo

for sample in $(ls data/process/metatrans_contigs/); do

        cp data/process/metatrans_contigs/"$sample"/transcripts.fasta "$scaffold_path"/"$sample"_scaffolds.fasta
        sed 's/^>/>'"$sample"'_/g' "$scaffold_path"/"$sample"_scaffolds.fasta >> "$scaffold_path"/all_raw_scaffolds.fasta

done


### Remove contigs under 1000 bases

echo "Removing contigs under 1000 bases"
echo

reformat.sh minlength=1000 overwrite=t in="$scaffold_path"/all_raw_scaffolds.fasta out="$scaffold_path"/all_long_scaffolds.fasta


### Create temporary contig file and remove line breaks from the sequences

echo "Creating temporary contig file and removing line breaks from the sequences"
echo
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' "$scaffold_path"/all_long_scaffolds.fasta > "$scaffold_path"/all_scaffolds.fasta

### Search for Circular Contigs

julia ccontigs.jl -i "$scaffold_path"/all_scaffolds.fasta -o "$cir_contig_path"/ccontigs_vlp_out.tsv

### Build Bowtie index files

bowtie2-build "$scaffold_path"/all_scaffolds.fasta data/references/metatrans_index/all_scaffolds 

### Run Read mappings, convert to BAM format, and sort

for line in $(cat data/process/trimmed_metatrans/samples.txt); do

	bowtie2 -x data/references/metatrans_index/all_scaffolds -1 data/process/trimmed_metatrans/"$line"_forward_paired.fastq\
	-2 data/process/trimmed_metatrans/"$line"_reverse_paired.fastq -U data/process/trimmed_metatrans/"$line"_unpaired.fastq\
	-S data/process/mapping/metatrans/"$line".sam --end-to-end --sensitive -p 12 -t > data/process/mapping/metatrans/"$line"_mapping_stats.txt

	samtools view -S -b data/process/mapping/metatrans/"$line".sam > data/process/mapping/metatrans/"$line".bam

	samtools sort data/process/mapping/metatrans/"$line".bam -o data/process/mapping/metatrans/"$line"_sorted.bam

	samtools idxstats data/process/mapping/metatrans/"$line"_sorted.bam > data/process/mapping/metatrans/"$line".tsv

done

