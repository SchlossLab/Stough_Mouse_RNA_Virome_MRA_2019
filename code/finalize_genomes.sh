#!/usr/bin/bash

grep -E "RNA-directed RNA polymerase|RNA-dependent RNA polymerase" data/raw/genome_candidates/annotations/genome_cand_annotations.tsv | awk '{ print $1 }' | uniq
