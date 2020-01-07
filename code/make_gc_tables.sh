#!/usr/bin/bash

for group in astrovirus mitovirus; do

	python code/calc_gc.py data/process/final_genomes/$group/fasta/ results/tables/"$group"_gc.csv

done


