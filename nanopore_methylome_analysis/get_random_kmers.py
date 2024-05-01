########## packages 
import Bio
from Bio import SeqIO
import random
import argparse
from argparse import RawTextHelpFormatter

##### Arguments #####

parser = argparse.ArgumentParser(
    description='input a ref sequence, specicy length of kmer and number of seqs you want. Returns randomly extracted kmers in multifasta format'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-ref', default = None, help='reference genome in fasta or multi fasta format, if multi, first seq will be used to extract kmers')
parser.add_argument('-k', type=int, help='length of kmer')
parser.add_argument('-n', type=int, help='number of random kmers to generate')
args = parser.parse_args()

#import reference genome and make all NTs uppercase
ref = list(SeqIO.parse(args.ref, "fasta"))
for i in range(len(ref)):
  ref[i].seq = ref[i].seq.upper()

#create multi fasta for control seeks
window = args.k #length of control seek
z = 1
for i in range(args.n): 
  position = random.randint(1, len(ref[0].seq)-11)
  seq = ref[0].seq[position:position+window]
  print('>seq_', z,'\n', seq, sep='')
  z += 1
