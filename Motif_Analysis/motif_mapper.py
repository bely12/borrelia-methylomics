#!/usr/bin/env python

import pandas as pd
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from tools_and_utilities import seq_tools
import argparse
from argparse import RawTextHelpFormatter

#ARGUMENTS
parser = argparse.ArgumentParser(
    description='Takes a list of motifs (1 per line, compatable with ambiguous NTs, and maps them to a genome).\n'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-fasta', help='fasta file of sequences')
parser.add_argument('-motif', help='file containing 1 motif per line')
parser.add_argument('-out', help='pre-fix for output file', default=None)
args = parser.parse_args()

#define the motifs to search for
with open(args.motif, 'r') as motifs:
  motif_list = motifs.readlines()
for i in range(len(motif_list)):
  motif_list[i] = motif_list[i].replace("\n", "")

#map motifs
mapped_motifs = seq_tools.motif_mapper(motif_list, args.fasta)

#print results
for item in mapped_motifs:
  print(item['id'],'\t',item['motif'],'\t',item['seq'],'\t',item['strand'],'\t',item['pos'])

# df = pd.DataFrame(mapped_motifs)
# df.to_csv(args.out+'_mapped_motifs.tsv', sep='\t', index=False, header=True)


