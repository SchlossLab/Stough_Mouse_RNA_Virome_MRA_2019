#!/usr/bin/bash
#SBATCH --job-name=blast_all_rdrp_10_04_2019
#SBATCH --mail-user=jmastough@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --account=pschloss
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --mem=48g

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
