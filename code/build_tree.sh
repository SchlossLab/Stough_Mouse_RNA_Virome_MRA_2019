#!/usr/bin/bash

####################################################################
# This script prepares raw short reads for assembly and mapping.   #
# To run this script you must have MAFFT, trimal, and iqtree       #
# installed.                                                       #
####################################################################

### Set up environment

tree_path="data/process/trees"
orf_path="data/process/final_genomes"
mito_ref_file="data/references/mito_rdrp_refs.fasta"
astro_ref_file="data/references/astrovirus_rdrp_refs.fasta"

### This chunk checks for whether the astrovirus tree directory exists
# and creates it if it doesn't. This directory is necessary for the
# trimming step in the code below.

for folder in trees fasta orfs; do 

	for taxa in astrovirus mitovirus; do

		if [ -d $tree_path/$taxa/$folder ]
		then
        		echo "$taxa tree folder already exists, continuing..."
        		echo
		else
        		echo "$taxa tree folder doesn't exist, creating and continuing..."
        		echo
        		mkdir -p $tree_path/$taxa/$folder
		fi

	done

done

### This chunk

for hit in $(grep -E "RNA-directed RNA polymerase" $orf_path/astrovirus/annotations/cef630_NODE_270_annotations.tsv | awk '{ print $1 }' | uniq); do

	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $orf_path/astrovirus/orfs/cef630_NODE_270_orf_translations.fasta > $orf_path/astrovirus/orfs/cef630_NODE_270_orf_translations_temp.fasta
        grep -A 1 "$hit" $orf_path/astrovirus/orfs/cef630_NODE_270_orf_translations_temp.fasta > $tree_path/astrovirus/orfs/$hit.fasta
	rm $orf_path/mitovirus/orfs/cef630_NODE_270_orf_translations_temp.fasta	

done

for genome in cefmock_NODE_3040 clinmock_NODE_4406 gf_NODE_1298 strep630_NODE_11363 strepmock_NODE_4960; do

	for hit in $(grep -E "RNA-dependent RNA polymerase" $orf_path/mitovirus/annotations/"$genome"_annotations.tsv | awk '{ print $1 }' | uniq); do
		
		awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $orf_path/mitovirus/orfs/"$genome"_orf_translations.fasta > $orf_path/mitovirus/orfs/"$genome"_orf_translations_temp.fasta
		grep -A 1 "$hit" $orf_path/mitovirus/orfs/"$genome"_orf_translations_temp.fasta > $tree_path/mitovirus/orfs/$hit.fasta
		rm $orf_path/mitovirus/orfs/"$genome"_orf_translations_temp.fasta
	done

done


### Combine contigs with refs

cat $mito_ref_file $tree_path/mitovirus/orfs/*.fasta > $tree_path/mitovirus/fasta/mito_rdrp_seqs_temp.fasta
sed -i 's/cefmock_NODE_3040_1/JS1/g' $tree_path/mitovirus/fasta/mito_rdrp_seqs_temp.fasta
sed -i 's/clinmock_NODE_4406_1/JS2/g' $tree_path/mitovirus/fasta/mito_rdrp_seqs_temp.fasta
sed -i 's/gf_NODE_1298_1/JS3/g' $tree_path/mitovirus/fasta/mito_rdrp_seqs_temp.fasta
sed -i 's/strep630_NODE_11363_1/JS4/g' $tree_path/mitovirus/fasta/mito_rdrp_seqs_temp.fasta
sed -i 's/strepmock_NODE_4960_1/JS5/g' $tree_path/mitovirus/fasta/mito_rdrp_seqs_temp.fasta
sed -i 's/NC_/NC/g' $tree_path/mitovirus/fasta/mito_rdrp_seqs_temp.fasta
cut -d ' ' -f 1 $tree_path/mitovirus/fasta/mito_rdrp_seqs_temp.fasta > $tree_path/mitovirus/fasta/mito_rdrp_seqs.fasta

cat $astro_ref_file $tree_path/astrovirus/orfs/*.fasta > $tree_path/astrovirus/fasta/astro_rdrp_seqs.fasta
sed -i 's/cef630_NODE_270_2/Murine astrovirus JS1/g' $tree_path/astrovirus/fasta/astro_rdrp_seqs.fasta
sed -i 's/ /_/g' $tree_path/astrovirus/fasta/astro_rdrp_seqs.fasta

### Alignment

mafft --globalpair --maxiterate 1000 $tree_path/mitovirus/fasta/mito_rdrp_seqs.fasta > $tree_path/mitovirus/fasta/mito_rdrp_aligned.fasta

mafft --globalpair --maxiterate 1000 $tree_path/astrovirus/fasta/astro_rdrp_seqs.fasta > $tree_path/astrovirus/fasta/astro_rdrp_aligned.fasta

### Alignment trimming

trimal -in $tree_path/mitovirus/fasta/mito_rdrp_aligned.fasta -out $tree_path/mitovirus/fasta/mito_rdrp_aligned_nogaps.fasta -gt 0.6

trimal -in $tree_path/astrovirus/fasta/astro_rdrp_aligned.fasta -out $tree_path/astrovirus/fasta/astro_rdrp_aligned_nogaps.fasta -gt 0.6

### Make Tree

iqtree -s $tree_path/mitovirus/fasta/mito_rdrp_aligned_nogaps.fasta -pre mito_tree -st AA -m MFP -bb 1000 -nt 12
mv mito_tree* $tree_path/mitovirus/trees/

iqtree -s $tree_path/astrovirus/fasta/astro_rdrp_aligned_nogaps.fasta -pre astro_tree -st AA -m MFP -bb 1000 -nt 12
mv astro_tree* $tree_path/astrovirus/trees/

