#!/usr/bin/bash

### PBS Preamble begin

#PBS -N trim_metatranscriptomes_5_29_2019
#PBS -l nodes=1:ppn=12
#PBS -l pmem=1gb
#PBS -l walltime=48:00:00
#PBS -A pschloss_fluxod
#PBS -q fluxod
#PBS -M jmastough@gmail.com
#PBS -m abe
#PBS -j oe
#PBS -V

### PBS Preamble End

### Change to the directory you submitted from
if [ -n "$PBS_O_WORKDIR" ]; then
        cd $PBS_O_WORKDIR
else
    	PBS_O_WORKDIR=$(pwd)
fi
echo "Job working directory:"
pwd
echo

source activate mj_omics

### Assemble trimmed reads

for sample in $(cat ./data/process/raw_unzipped_metatrans/samples.txt); do

        trimmomatic PE -threads 12 data/process/raw_unzipped_metatrans/"$sample"_forward.fastq data/process/raw_unzipped_metatrans/"$sample"_reverse.fastq data/process/trimmed_metatrans/"$sample"_forward_paired.fastq data/process/trimmed_metatrans/"$sample"_forward_unpaired.fastq data/process/trimmed_metatrans/"$sample"_reverse_paired.fastq data/process/trimmed_metatrans/"$sample"_reverse_unpaired.fastq LEADING:10 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:50 ILLUMINACLIP:data/references/TruSeq3-PE.fa:2:30:10
	cat data/process/trimmed_metatrans/"$sample"_forward_unpaired.fastq data/process/trimmed_metatrans/"$sample"_reverse_unpaired.fastq > data/process/trimmed_metatrans/"$sample"_unpaired.fastq

done


