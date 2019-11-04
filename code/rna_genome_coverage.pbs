#! /usr/bin/bash

### PBS Preamble begin

#PBS -N rna_genome_coverage_10_25_2019
#PBS -l nodes=1:ppn=12
#PBS -l pmem=4gb
#PBS -l walltime=12:00:00
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


### Build genome histograms

for file in $(ls data/process/metatrans_contigs); do

	bedtools coverage -a data/process/mapping/metatrans/genome_tables/"$file"_blasts.tsv -b data/process/mapping/metatrans/"$file"_sorted.bam > data/process/mapping/metatrans/coverage/"$file"_coverage.tsv

done
