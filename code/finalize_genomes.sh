#!/usr/bin/bash

gen_dir="data/process/final_genomes"
cand_dir="data/raw/genome_candidates"

### This chunk checks for whether the finalized genome directories 
# exist and creates them if they don't. These directories are 
# necessary for the steps in the code below.

if [ -d $gen_dir ]
then
    	echo "Genomes folder already exists, continuing..."
        echo
else
    	echo "Genomes folder doesn't exist, creating and continuing..."
        echo
	mkdir $gen_dir
fi

for folder in mitovirus astrovirus; do

	if [ -d "$gen_dir"/"$folder" ]
	then
    		echo "$folder already exists, continuing..."
        	echo
	else
    		echo "$folder folder doesn't exist, creating and continuing..."
        	echo
		mkdir "$gen_dir"/"$folder"     
	fi

done

for folder in $(ls $gen_dir); do

	for subdir in fasta orfs annotations genbank; do

		if [ -d "$gen_dir"/"$folder"/"$subdir" ]
		then
    			echo "$folder $subdir already exists, continuing..."
        		echo
		else
    			echo "$folder $subdir folder doesn't exist, creating and continuing..."
        		echo
			mkdir "$gen_dir"/"$folder"/"$subdir"     
		fi
	
	done

done


### This chunk 

for hit in $(grep -E "RNA-directed RNA polymerase" $cand_dir/annotations/genome_cand_annotations.tsv | awk '{ print $1 }' | uniq | cut -d '_' -f 1,2,3); do

	grep -A 1 "$hit$" $cand_dir/genomes/rdrp_hits.fasta > $gen_dir/astrovirus/fasta/$hit.fasta

done

for hit in $(grep -E "RNA-dependent RNA polymerase" $cand_dir/annotations/genome_cand_annotations.tsv | awk '{ print $1 }' | uniq | cut -d '_' -f 1,2,3); do

        grep -A 1 "$hit$" $cand_dir/genomes/rdrp_hits.fasta > $gen_dir/mitovirus/fasta/$hit.fasta

done


### This chunk

for group in astrovirus mitovirus; do

	for genome in $(ls "$gen_dir"/"$group"/fasta/ | cut -d '.' -f 1); do

		prodigal -a $gen_dir/"$group"/orfs/"$genome"_orf_translations.fasta -c -d $gen_dir/"$group"/orfs/"$genome"_orf_nucleotides.fasta -f sco -g 1 -i "$gen_dir"/"$group"/fasta/"$genome".fasta -m -o "$gen_dir"/"$group"/orfs/"$genome"_orfs.sco -p meta
		sed -i 's/*//g' "$gen_dir/"$group"/orfs/$genome"_orf_translations.fasta
		interproscan.sh -dp -b "$gen_dir"/"$group"/annotations/"$genome"_annotations -cpu 12 -f tsv -i "$gen_dir"/"$group"/orfs/"$genome"_orf_translations.fasta -t p -ms 50

	done

done
