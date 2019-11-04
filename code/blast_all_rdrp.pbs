#! /usr/bin/bash

### PBS Preamble begin

#PBS -N blast_all_rdrp_10_04_2019
#PBS -l nodes=1:ppn=12
#PBS -l mem=48gb
#PBS -l walltime=36:00:00
#PBS -A pschloss_fluxod
#PBS -q fluxod
#PBS -M andbeaud@umich.edu
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
ref_path="data"

### Build marker gene databases

# Bacteriocin database
makeblastdb -in "$ref_path"/references/blast_refs/all_rdrp.fasta -dbtype prot -title all_rdrp -out "$ref_path"/references/blast_refs/all_rdrp

### Define a function that blasts orf translations against marker gene databases

for treatment in $(ls data/process/metatrans_contigs); do

	blastx -query data/process/scaffolds/metatrans_scaffolds/"$treatment"_transcripts.fasta -query_gencode 11 -db "$ref_path"/references/blast_refs/all_rdrp -evalue 1e-20 -out data/process/blasts/all_rdrp/"$treatment"_blasts.tsv -outfmt 6 

done
