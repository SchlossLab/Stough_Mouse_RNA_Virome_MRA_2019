#!/usr/bin/bash

############################################################################################
# This script extracts RdRP hit scaffolds using the blast tables produced during           #
# blast_allrdrp.sh and deposits them into the data/process/genome_candidates/genomes       #
# directory. The script then uses prodigal to predict open reading frames and deposits     #
# the results ass gff files in the data/process/genome_candidates/orfs directory. Prodigal #
# also saves the protein translations in fasta format to the same directory. The script    #
# then deletes * characters from the translations and runs interproscan to annotate        #
# protein function, depositing the results in the                                          #
# data/process/genome_candidates/annotations directory.                                    # 
# This script requires a working installation of the Prodigal Open Reading Frame           #
# prediction software.                                                                     #
###########################################################################################

### This chunk checks for whether the genome_candidates data directory exists
# and creates it if it doesn't. This directory is necessary for the
# trimming step in the code below.

if [ -d "data/raw/genome_candidates/" ]
then
    	echo "Genome_candidates folder already exists, continuing..."
        echo
else
    	echo "Genome_candidates folder doesn't exist, creating and continuing..."
        echo
	mkdir data/raw/genome_candidates
fi

### This chunk checks for whether the genomes data directory exists
# and creates it if it doesn't. This directory is necessary for the
# trimming step in the code below.

if [ -d "data/raw/genome_candidates/genomes" ]
then
    	echo "Genomes folder already exists, continuing..."
        echo
else
    	echo "Genomes folder doesn't exist, creating and continuing..."
        echo
	mkdir data/raw/genome_candidates/genomes
fi

### This chunk checks for whether the orfs data directory exists
# and creates it if it doesn't. This directory is necessary for the
# trimming step in the code below.

if [ -d "data/raw/genome_candidates/orfs" ]
then
        echo "ORFs folder already exists, continuing..."
        echo
else
    	echo "ORFs folder doesn't exist, creating and continuing..."
        echo
        mkdir data/raw/genome_candidates/orfs
fi

### This chunk checks for whether the annotations data directory exists
# and creates it if it doesn't. This directory is necessary for the
# trimming step in the code below.

if [ -d "data/raw/genome_candidates/annotations" ]
then
    	echo "Annotations folder already exists, continuing..."
        echo
else
    	echo "Annotations folder doesn't exist, creating and continuing..."
        echo                                                                   
    	mkdir data/raw/genome_candidates/annotations
fi

### Extract RNA virus genomes from contig files

if [ -e data/raw/genome_candidates/genomes/rdrp_hits.fasta ];
then
    	echo "Deleting old hits file"
        rm data/raw/genome_candidates/genomes/rdrp_hits.fasta
else
    	echo "No old contig files to delete"
fi

for hit in $(awk '{ print $1 }' data/raw/blasts/rdrp_blasts.tsv | uniq); do

	grep -A 1 "$hit" data/raw/scaffolds/all_scaffolds.fasta >> data/raw/genome_candidates/genomes/rdrp_hits.fasta

done

grep ">" data/raw/genome_candidates/genomes/rdrp_hits.fasta | sed 's/>//g' > results/tables/rdrp_hits.tsv

### This chunk uses Prodigal to predict open reading frames in the scaffolds with hits
# in the RdRP protein database and outputs nucleotide sequences, protein translations,
# and ORF locations in GFF format in the data/raw/genome_candidates/orfs/ directory.

prodigal -a data/raw/genome_candidates/orfs/rdrp_hit_orf_translations.fasta -c -d data/raw/genome_candidates/orfs/rdrp_hit_orf_nucleotides.fasta -f gff -g 1 -i data/raw/genome_candidates/genomes/rdrp_hits.fasta -m -o data/raw/genome_candidates/orfs/rdrp_hit_orfs.gff -p meta


### This chunk removes asterisks from the ORF translation output from prodigal. 
# Asterisks represent stop codons and interproscan will fail if they're not
# removed.

sed -i 's/*//g' data/raw/genome_candidates/orfs/rdrp_hit_orf_translations.fasta


### This chunk uses interproscan to annotate ORFs output by prodigal and output the
# results in gff3 and tsv format in the data/raw/genome_candidates/annotations/
# directory.

interproscan.sh -dp -b data/raw/genome_candidates/annotations/genome_cand_annotations -cpu 12 -f tsv,gff3 -i data/raw/genome_candidates/orfs/rdrp_hit_orf_translations.fasta -t p -ms 50

