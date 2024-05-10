#!/usr/bin/env python

import seq_tools
import argparse
from argparse import RawTextHelpFormatter

########## ARGUMENTS ##########
parser = argparse.ArgumentParser(
    description='Does what the name says it does.\n'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-fasta', help='path to your fasta file\n'
                    ' \n')

parser.add_argument('-n', type=int, default=1, help='how many different times you want to shuffle your sequence(s)\n'
                    ' \n')

args = parser.parse_args()

sim_seqs = seq_tools.syn_codon_shuffle(args.fasta, args.n)

for item in sim_seqs:
  print('>',item['id'],'_sim',item['num'],'\n',item['seq'],sep='')