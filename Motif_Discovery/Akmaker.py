########## packages 
import Bio
from Bio import SeqIO
import pandas as pd
import random
from random import sample
import math
import seq_tools
import argparse
from argparse import RawTextHelpFormatter

##### Arguments #####

parser = argparse.ArgumentParser(
    description=''
    ' \n', formatter_class=RawTextHelpFormatter)


### Arguments
parser.add_argument('-bed', help='tsv bed file with mod calls and kmers. kmers must be in the 3rd column of the table, no header\n')
parser.add_argument('-ref', default = None, help='reference genome in fasta or multi fasta format, first seq will be used to generate control seqs')

parser.add_argument('-all', const=True, nargs='?')
parser.add_argument('-sample', type=int, help='number of random adenine kmers if selecting all arg', default=None)
parser.add_argument('-modified', const=True, nargs='?')

parser.add_argument('-len', default=11, type=int, help='size of kmer you wish to extract')
parser.add_argument('-depth', default=10, type=int, help='minimum sequencing depth to use, inclusive')
parser.add_argument('-mod_threshold', default=50.0, type=float, help='percentage of reads that need to be predicted as modified for that position to call position 6mA')

parser.add_argument('-controls', const=True, nargs='?', help='outputs a set of control sequences that are generated ')
parser.add_argument('-out', help='prefix for output fasta file')

parser.add_argument('-random', const=True, nargs='?')
parser.add_argument('-n', type=int, help='number of random kmers to generate')
args = parser.parse_args()


# import reference genome and make all NTs uppercase
ref = list(SeqIO.parse(args.ref, "fasta"))
for i in range(len(ref)):
  ref[i].seq = ref[i].seq.upper()

if args.bed != None:
  # import bedmethyl table and do some formatting
  bed = pd.read_table(args.bed, sep = '\t', header=None)
  bed[['cov', 'percent_mod', 'N_mod', 'N_can', 'N_other', 'N_del', 'N_fail', 'N_diff', 'N_nocall']] = bed[9].str.split(' ', expand=True)
  bed = bed.rename(columns={0: 'chromosome', 1: 'start', 2: 'end', 3: 'mod', 4: 'score', 5: 'strand'})
  bed = bed.drop(columns=[6,7,8,9])
  bed['percent_mod'] = bed['percent_mod'].astype(float)
  bed['cov'] = bed['cov'].astype(int)

  # random sample adenine kmers
  if args.all == True:
    df = bed[ bed['cov'] >= args.depth] #filter for positions with at least 10X depth 
    recs = df.to_dict('records')
    if args.sample != None:
      recs = sample(recs, args.sample)

  # get kmers for likely 6mA positions
  if args.modified == True:
    df = bed[ (bed['percent_mod'] > args.mod_threshold) & (bed['cov'] >= args.depth)]
    recs = df.to_dict('records')

  # extract 11mers for every adenine 
  kmers = []
  down = math.floor(args.len/2)
  up = math.ceil(args.len/2)
  for i in range(len(recs)):
    plasmid = recs[i]['chromosome']
    position = recs[i]['start']
    for j in range(len(ref)):
      if ref[j].id == plasmid:
        if (position < 100) or (position > len(ref[j].seq)-100): #to not select from extreme ends of seqs 
          break
        #seq = str(ref[j].seq[position-5 : position+6]) #extract 11mer with A in middle
        seq = str(ref[j].seq[position-down : position+up]) #for any size kmer
        if recs[i]['strand'] == '-':
          seq = seq_tools.reverse_complement(seq)
        kmers.append(seq)
        break

  # remove errant kmers that do not have A in middle position
  #fixed_kmers = [i for i in kmers if not (i[5] != 'A')]
  fixed_kmers = [i for i in kmers if not (i[down] != 'A')] #for any size kmer

  # print results in fasta format
  k = 1
  for seq in fixed_kmers:
    if args.all != None:
      print('>seq_',k,'\n',seq, sep='', file = open(args.out+'_sampled_adenine_seqs.fasta', "a"))
      k += 1
    else:
      print('>seq_',k,'\n',seq, sep='', file = open(args.out+'_positive_seqs.fasta', "a"))
      k += 1

  if args.controls == True:
    #create multi fasta for control seeks
    z = 1
    for i in range(k): 
      position = random.randint(1, len(ref[0].seq)-args.len)
      seq = ref[0].seq[position:position+args.len]
      print('>control_', z,'\n', seq, sep='', file = open(args.out+'_control_seqs.fasta', "a"))
      z += 1

# grab random kmers from ref genome and make fasta
if args.random != None:
  window = args.len #length of control seek
  z = 1
  for i in range(args.n): 
    position = random.randint(1, len(ref[0].seq)-args.len)
    seq = ref[0].seq[position:position+window]
    print('>seq_', z,'\n', seq, sep='')
    z += 1
