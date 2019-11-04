#!/usr/bin/bash

### PBS Preamble begin

#PBS -N predict_orfs_10_25_2019
#PBS -l nodes=1:ppn=12
#PBS -l pmem=4gb
#PBS -l walltime=06:00:00
#PBS -A pschloss_fluxod
#PBS -q fluxod
#PBS -M jmastough@gmail.com
#PBS -m abe
#PBS -j oe
#PBS -V

### Change to the directory you submitted from
if [ -n "$PBS_O_WORKDIR" ]; then
        cd $PBS_O_WORKDIR
else
    	PBS_O_WORKDIR=$(pwd)
fi
echo "Job working directory:"
pwd
echo

### Set up environment

source activate mj_omics

### Extract RNA virus genomes from contig files

for table in $(ls data/process/blasts/all_rdrp/); do

	awk '{ print $1 }' data/process/blasts/all_rdrp/"$table" > data/process/blasts/rdrp_hits/"$table"

done

cat data/process/blasts/rdrp_hits/* | uniq > data/process/blasts/rdrp_hits/all_hits.tsv

if [ -e data/process/rna_virus_genomes/genomes/all_rdrp_hits.fasta ];
then
    	echo "Deleting old hits file"
        rm data/process/rna_virus_genomes/genomes/all_rdrp_hits.fasta
else
    	echo "No old contig files to delete"
fi


for id in $(cat data/process/blasts/rdrp_hits/all_hits.tsv); do

	grep -A 1 "$id" data/process/scaffolds/metatrans/all_temp_scaffolds.fasta >> data/process/rna_virus_genomes/genomes/all_rdrp_hits.fasta

done

### Predict ORFS

prodigal -a data/process/rna_virus_genomes/orfs/rdrp_hit_orf_translations.fasta -c -d data/process/rna_virus_genomes/orfs/rdrp_hit_orf_nucleotides.fasta -f gff -g 1 -i data/process/rna_virus_genomes/genomes/all_rdrp_hits.fasta -m -o data/process/rna_virus_genomes/orfs/rdrp_hit_orfs.gff -p meta
