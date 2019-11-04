#! /usr/bin/bash

### PBS Preamble begin

#PBS -N assemble_metatrans_6_26_2019
#PBS -l nodes=1:largemem:ppn=40
#PBS -l mem=900gb
#PBS -l walltime=500:00:00
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

### Create temporary directory, which is removed at end of run
TEMP=$(mktemp -d)
if [ $? -ne 0 ]
then
    	echo "can't create $TEMP ... exit"
        exit 1
fi

### Set up environment

source activate mj_omics
read_path="$HOME/mj_omics/data/process/trimmed_metatrans"
contig_path="$HOME/mj_omics/data/process/metatrans_contigs"

### Assemble trimmed reads

for file in $(cat ./data/process/trimmed_metatrans/samples.txt); do
	mkdir "$contig_path"/"$file"
	spades.py --rna -o "$contig_path"/"$file"/ -1 "$read_path"/"$file"_forward_paired.fastq -2 "$read_path"/"$file"_reverse_paired.fastq -s "$read_path"/"$file"_unpaired.fastq -t 40 -m 900 --tmp-dir $TEMP
done

### Clean up temp folder
rm -rf $TEMP
