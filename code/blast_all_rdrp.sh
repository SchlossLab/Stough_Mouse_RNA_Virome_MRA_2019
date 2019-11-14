#!/usr/bin/bash

###########################################################################
# This script uses blastx to screen scaffolds produced during assembly    #
# against a database of RefSeq RNA-dependent RNA Polymerase protein       #
# sequences.                                                              #
# A working installation of NCBI BLAST+ is required to to run this script.#
###########################################################################

### This chunk checks for whether the blasts directory exists
# and creates it if it doesn't. This directory is necessary for the
# steps in the code below.

if [ -d "data/process/blasts/" ]
then
    	echo "Blasts folder already exists, continuing..."
        echo
else
    	echo "Blasts folder doesn't exist, creating and continuing..."
        echo
	mkdir data/process/blasts
fi


### This chunk uses makeblastdb from the BLAST+ software suite to build
# a blast database from the all_rdrp.fasta file containing RefSeq
# RdRP protein sequences and deposits it in the data/references directory

makeblastdb -in data/references/all_rdrp.fasta -dbtype prot -title all_rdrp -out data/references/all_rdrp


### This chunk blasts processed scaffolds against the RdRP database
# produced above and deposts the results in tsv format in the 
# data/process/blasts/ directory.

blastx -query data/raw/scaffolds/all_scaffolds.fasta -query_gencode 11 -db data/references/all_rdrp -evalue 1e-20 -out data/process/blasts/rdrp_blasts.tsv -outfmt 6 

