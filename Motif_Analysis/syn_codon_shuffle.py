#!/usr/bin/env python

import seq_tools
import argparse
from argparse import RawTextHelpFormatter

#ARGUMENTS
parser = argparse.ArgumentParser(
    description='Does what the name says it does. Save on command line using > filename.txt\n'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-fasta', help='path to your fasta file\n')
parser.add_argument('-n', type=int, default=1, help='how many different times you want to shuffle your sequence(s)\n')
args = parser.parse_args()

#synonymous codon shuffle
sim_seqs = seq_tools.syn_codon_shuffle(args.fasta, args.n)

#print output in fasta format
for item in sim_seqs:
  print('>',item['id'],'_sim',item['num'],'\n',item['seq'],sep='')