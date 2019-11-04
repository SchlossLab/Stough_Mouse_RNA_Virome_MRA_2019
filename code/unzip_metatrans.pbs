#! /usr/bin/bash

### PBS Preamble begin

#PBS -N cat_raw_reads_6_5_2019
#PBS -l nodes=1:ppn=12
#PBS -l pmem=1gb
#PBS -l walltime=12:00:00
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

### Unzip fastq files

for folder in $(ls ./data/raw/metatrans_raw/); do

	gunzip ./data/raw/metatrans_raw/"$folder"/*.fastq.gz
	cat ./data/raw/metatrans_raw/"$folder"/*R1_0??.fastq > ./data/process/raw_unzipped_metatrans/"$folder"_forward.fastq
	cat ./data/raw/metatrans_raw/"$folder"/*R2_0??.fastq > ./data/process/raw_unzipped_metatrans/"$folder"_reverse.fastq
	
done
