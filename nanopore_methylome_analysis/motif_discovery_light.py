### packages ###
import Bio
from Bio import SeqIO
import seq_tools
import mod_tools
import argparse
from argparse import RawTextHelpFormatter

### ARGUMENTS ###
parser = argparse.ArgumentParser(
    description='takes an input of short candidate seqs with adenine in center, outputs n most common kmers of spec length'
    ' \n', formatter_class=RawTextHelpFormatter)

def list_of_ints(arg):
    return list(map(int, arg.split(',')))

parser.add_argument('-candidates', help='mutli fasta file with candidate kmers centered on mod adenine')
parser.add_argument('-mod_pos_index', help='index position of mod adenine in candidate seqs', type=int)
parser.add_argument('-target_lengths', help='list containing the lengths of kmers to test to find targets', type=list_of_ints)
parser.add_argument('-n_targets', type=int, help='number of potential targets to record based on top counts for each target length tested', default=15)
args = parser.parse_args()


#load candidate motifs (output from get_adenine_kmers.py)
seqList = list(SeqIO.parse(args.candidates, "fasta"))

#find potential target motifs
test_motifs = mod_tools.top_kmers(seqList, args.mod_pos_index, args.target_lengths, args.n_targets)

#print results
for sets in test_motifs:
  for i in range(len(sets)):
    print(sets[i]['kmer_len'],'\t',sets[i]['seq'],'\t',sets[i]['count'],'\t',sets[i]['mod_pos'])