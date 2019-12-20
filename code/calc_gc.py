#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 13:24:27 2019

@author: stoughj
"""

from Bio import SeqIO
import pandas as pd
import os
import sys

input_dir = sys.argv[1]
input_file = os.listdir(input_dir)
input_seqs = [SeqIO.read(f"{input_dir}{file}", "fasta") for file in input_file]

def calc_gc(input_seq_obj):
    g_count = input_seq_obj.seq.count("G")
    c_count = input_seq_obj.seq.count("C")
    gc = ((g_count + c_count)/len(input_seq_obj.seq))*100
    return gc

id_gc = {genome.id:calc_gc(genome) for genome in input_seqs}
output_table = pd.DataFrame(id_gc.items(), columns = ["contig_id", "gc_content"])

output_table.to_csv(sys.argv[2], index = False)

#