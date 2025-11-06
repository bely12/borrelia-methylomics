########## packages 
import Bio
from Bio import SeqIO
import pandas as pd
import numpy as np
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from tools_and_utilities import seq_tools
from tools_and_utilities import mod_tools
import argparse
from argparse import RawTextHelpFormatter

### ARGUMENTS ###
parser = argparse.ArgumentParser(
    description='calculates what percentage of a specified motif is modified in the genome'
    ' \n', formatter_class=RawTextHelpFormatter)

def list_of_ints(arg):
    return list(map(int, arg.split(',')))

parser.add_argument('-bed', help='tsv bed file with mod calls and kmers. kmers must be in the 3rd column of the table, no header\n')
parser.add_argument('-ref', default = None, help='reference genome in fasta or multi fasta format, first seq will be used to generate control seqs')
parser.add_argument('-candidates', help='mutli fasta file with candidate kmers centered on mod adenine')
parser.add_argument('-mod_base', default = 'A', help='what modified base ex. A or C')
parser.add_argument('-mod_pos_index', help='index position of mod adenine in candidate seqs', type=int)
parser.add_argument('-mod_call_thresh', type=float, default=50.0, help='threshold for calling a position modified in bed file')
parser.add_argument('-min_cov', type=int, default=10, help='mininum position coverage to filter for')
parser.add_argument('-target_lengths', help='list containing the lengths of kmers to test to find targets', type=list_of_ints)
parser.add_argument('-n_targets', type=int, help='number of potential targets to record based on top counts for each target length tested', default=15)
args = parser.parse_args()


#import reference genome and change all letters to uppercase
ref = list(SeqIO.parse(args.ref, "fasta"))
for i in range(len(ref)):
  ref[i].seq = ref[i].seq.upper()

#prepare bed file and convert to list of row dictionaries 
bed = pd.read_table(args.bed, sep = '\t', header=None)
#bed[['cov', 'mod_freq', 'N_mod', 'N_can', 'N_other', 'N_del', 'N_fail', 'N_diff', 'N_nocall']] = bed[9].str.split(' ', expand=True)
#bed = bed.rename(columns={0: 'chromosome', 1: 'start', 2: 'end', 3: 'mod', 4: 'score', 5: 'strand'})
#bed = bed.drop(columns=[6,7,8,9])
bed = bed.rename(columns={0: 'chromosome', 
                          1: 'start', 
                          2: 'end', 
                          3: 'mod', 
                          4: 'score', 
                          5: 'strand',
                          9:'cov',
                          10:'mod_freq',
                          11:'N_mod',
                          12:'N_can',
                          13:'N_other',
                          14:'N_del',
                          15:'N_fail',
                          16:'N_diff',
                          17:'N_nocall'})
bed = bed.drop(columns=[6,7,8])
bed['mod_freq'] = bed['mod_freq'].astype(float)
bed['cov'] = bed['cov'].astype(int)
recs = bed.to_dict('records')

#get 11mers for each position
bed2 = []
#for every row in bedmethyl table
for i in range(len(recs)):
  plasmid = recs[i]['chromosome'] #specify plasmid
  position = recs[i]['start'] #get coordinate of adenine
  #find the plasmid in the reference genome and extraxt kmer
  for j in range(len(ref)):
    if ref[j].id == plasmid:
      #seq = str(ref[j].seq[position-5 : position+6]) # for 11mers
      seq = str(ref[j].seq[position-15 : position+16]) # for longer 31mers
      if recs[i]['strand'] == '-':
        seq = seq_tools.reverse_complement(seq)
      #make a new list of dictionaries by adding 11mer as key value pair to existing bedmethyl dictionaries
      bed2.append(dict(recs[i], kmer = seq))
      break
#bed2

#remove errant kmers (don't have potential 6mA in correct position)
#bed3 = [i for i in bed2 if not ( (len(i['kmer']) != 11) or (i['kmer'][5] != 'A') )] # for 11mers
bed3 = [i for i in bed2 if not ( (len(i['kmer']) != 31) or (i['kmer'][15] != args.mod_base) )] # for 31mers

### get most common kmers in candidate seqs ###
#load candidate motifs (output from get_adenine_kmers.py)
seqList = list(SeqIO.parse(args.candidates, "fasta"))
#find potential target motifs
test_motifs = mod_tools.top_kmers(seqList, args.mod_pos_index, args.target_lengths, args.n_targets, mod_base = args.mod_base)


### caclculate percent modified for the most common kmers within candidate seqs ###
#set column name for results table 
print('motif','\t','occurences','\t','mod','\t','percent')

#calculate
for sets in test_motifs:
  for i in range(len(sets)):
    motif_list = [sets[i]['seq']]
    mod_pos = sets[i]['mod_pos']
    results = mod_tools.get_mod_frequency(motif_list, bed3, mod_pos = mod_pos, mod_call_thresh = args.mod_call_thresh, min_cov = args.min_cov)
    for i in range(len(results)):
      if results[i]['percent_mod'] > 0.09:
        print(results[i]['motif'],'\t',
              results[i]['occurences'],'\t',
              results[i]['modified'],'\t',
              results[i]['percent_mod'])