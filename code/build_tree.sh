#!/usr/bin/bash

####################################################################
# This script prepares raw short reads for assembly and mapping.   #
# To run this script you must have MAFFT, trimal, and iqtree       #
# installed.                                                       #
####################################################################

### Set up environment

tree_path = "data/process/narnaviridae_tree"
contig_file = "data/raw/genome_candidates/genomes/narna_rdrp_hits.fasta"
mito_ref_file = "data/references/mito_rdrp_refs.fasta"

### This chunk checks for whether the trimmed data directory exists
# and creates it if it doesn't. This directory is necessary for the
# trimming step in the code below.

if [ -d $tree_path ]
then
        echo "Narna tree folder already exists, continuing..."
        echo
else
        echo "Narna tree folder doesn't exist, creating and continuing..."
        echo
        mkdir $tree_path
fi

### Combine contigs with refs

cat $mito_ref_file $contig_file > $tree_path/mito_rdrp_seqs.fasta

### Alignment

mafft --globalpair --maxiterate 1000 $tree_path/mito_rdrp_seqs.fasta > $tree_path/mito_rdrp_aligned.fasta


### Alignment trimming

trimal -in $tree_path/mito_rdrp_aligned.fasta -out $tree_path/mito_rdrp_aligned_nogaps.fasta -gt 0.5


### Make Tree

iqtree -s $tree_path/mito_rdrp_aligned_nogaps.fasta -pre mito_tree -st AA -m MFP -bb 1000 -nt 12

### Edit tree




### Output tree figure
