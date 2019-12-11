#!usr/bin/bash

###########################################################################
# This script processes assembled scaffolds, adding treatment groups to   #
# fasta sequence headers, removing short contigs, removing duplicate      #
# contigs, removing line breaks from fasta sequences, and shortening      #
# fasta headers in the final scaffold file.                               #
# Working installations of BBMap and CContigs is required to              #
# to run this script.                                                     #
###########################################################################


### This chunk checks for whether the scaffold directory exists
# and creates it if it doesn't. This directory is necessary for the
# steps in the code below.

if [ -d "data/raw/scaffolds/" ]
then
    	echo "Scaffolds read folder already exists, continuing..."
        echo
else
    	echo "Scaffolds read folder doesn't exist, creating and continuing..."
        echo
        mkdir data/raw/scaffolds
fi

### This chunk adds filenames to fasta headers, concatenates scaffolds into a 
# single file, and removes scaffolds under 1000 bases long

echo "Removing contigs under 1kb, adding treatment groups to fasta headers, and concatenating contigs into single file"
echo

for treatment in $(awk '{ print $2 }' data/process/samples.tsv); do

        reformat.sh minlength=1000 overwrite=t in=data/raw/assembly/"$treatment"/transcripts.fasta out=data/raw/scaffolds/"$treatment"_scaffolds.fasta
        sed 's/^>/>'"$treatment"'_/g' data/raw/scaffolds/"$treatment"_scaffolds.fasta >> data/raw/scaffolds/all_long_scaffolds.fasta

done

### This chunk uses BBMap function dedupe.sh to remove duplicate contigs 
# assembled in different treatment groups, then deletes the temporary file.

dedupe.sh in=data/raw/scaffolds/all_long_scaffolds.fasta out=data/raw/scaffolds/all_scaffolds_nodupes.fasta outd=data/raw/scaffolds/duplicates.fasta

### This chunk uses grep to extract contig names and stats into a tsv file
# that will be parsed in R in another script

grep ">" data/raw/scaffolds/all_long_scaffolds.fasta | sed 's/>//g' > results/tables/contig_stats_raw.tsv


### This chunk removes line breaks from the sequences for easy extraction
# later using grep

echo "Removing line breaks from the scaffold sequences and creating whole scaffold file"
echo
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' data/raw/scaffolds/all_scaffolds_nodupes.fasta > data/raw/scaffolds/all_temp_scaffolds.fasta

### This chunk shortens the contig names using sed to replace the treatment
# group names with shorter versions, and cut to remove all contig stats
# except the id number

sed -i 's/>cefoperazone_630_/>cef630_/g' data/raw/scaffolds/all_temp_scaffolds.fasta
sed -i 's/>cefoperazone_mock_/>cefmock_/g' data/raw/scaffolds/all_temp_scaffolds.fasta
sed -i 's/>clindamycin_630_/>clin630_/g' data/raw/scaffolds/all_temp_scaffolds.fasta
sed -i 's/>clindamycin_mock_/>clinmock_/g' data/raw/scaffolds/all_temp_scaffolds.fasta
sed -i 's/>streptomycin_630_/>strep630_/g' data/raw/scaffolds/all_temp_scaffolds.fasta
sed -i 's/>streptomycin_mock_/>strepmock_/g' data/raw/scaffolds/all_temp_scaffolds.fasta
sed -i 's/>germ_free_/>gf_/g' data/raw/scaffolds/all_temp_scaffolds.fasta
cut -d '_' -f 1,2,3 data/raw/scaffolds/all_temp_scaffolds.fasta > data/raw/scaffolds/all_scaffolds.fasta

rm data/raw/scaffolds/all_long_scaffolds.fasta data/raw/scaffolds/all_scaffolds_nodupes.fasta data/raw/scaffolds/duplicates.fasta data/raw/scaffolds/all_temp_scaffolds.fasta
