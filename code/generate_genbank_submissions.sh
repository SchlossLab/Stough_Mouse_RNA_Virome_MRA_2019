#!/usr/bin/bash

for group in $(ls data/process/final_genomes); do

	for virus in $(ls data/process/final_genomes/"$group"/fasta | cut -d '.' -f 1); do

		Rscript code/build_genome_tables.R data/process/final_genomes/"$group"/orfs/"$virus"_orfs.sco data/process/final_genomes/"$group"/annotations/"$virus"_annotations.tsv
		cp data/process/final_genomes/"$group"/fasta/"$virus".fasta data/process/final_genomes/"$group"/submission
		cp data/process/final_genomes/"$group"/annotations/"$virus"_annotation_table.tsv data/process/final_genomes/"$group"/submission

	done

done

