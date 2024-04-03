########## packages 
import seq_tools
import Bio
from Bio import SeqIO
import pandas as pd
import random
import argparse
from argparse import RawTextHelpFormatter

########## ARGUMENTS
parser = argparse.ArgumentParser(
    description='generate negative control motifs to use as input (motif_file arg) in motif_discovery.py script \n'
    'prints output onto console, so can save as txt file by using >'
    ' \n', formatter_class=RawTextHelpFormatter)

#input files
parser.add_argument('-bed', help='tsv bed file with mod calls and kmers. kmers must be in the 3rd column of the table, no header\n')
parser.add_argument('-ref', default = None, help='reference genome in fasta or multi fasta format')

#characteristics of control motifs
parser.add_argument('-len', type=int, help='length of motif')
parser.add_argument('-mod_position', type=int, help='position of m6A in motif')
parser.add_argument('-n', type=int, help='number of control motifs to generate')

#positive target/s to avoid 
parser.add_argument('-target', default=None, help='positive target motif to avoid in neg control sequences')
parser.add_argument('-target_list', default=None, nargs='*', help='list of positive target motifs to avoid in neg control motifs, each separated by comma')

args = parser.parse_args()

#import reference genome and make all NTs uppercase
ref = list(SeqIO.parse(args.ref, "fasta"))
for i in range(len(ref)):
  ref[i].seq = ref[i].seq.upper()

#import bedmethyl table, format, and turn rows in list of dictionaries 
bed = pd.read_table(args.bed, sep = '\t', header=None)
bed[['cov', 'percent_mod', 'N_mod', 'N_can', 'N_other', 'N_del', 'N_fail', 'N_diff', 'N_nocall']] = bed[9].str.split(' ', expand=True)
bed = bed.rename(columns={0: 'chromosome', 1: 'start', 2: 'end', 3: 'mod', 4: 'score', 5: 'strand'})
bed = bed.drop(columns=[6,7,8,9])
bed['percent_mod'] = bed['percent_mod'].astype(float)
bed['cov'] = bed['cov'].astype(int)
df = bed[ bed['cov'] >= 10] #filter for positions with at least 10X depth 
recs = df.to_dict('records')


#extract seqences surrounding every adenine position as a list called kmers (creates 11mers with A in middle)
kmers = []
for i in range(len(recs)):
  plasmid = recs[i]['chromosome']
  position = recs[i]['start']
  for j in range(len(ref)):
    if ref[j].id == plasmid:
      if (position < 100) or (position > len(ref[j].seq)-100): #to not select from extreme ends of seqs 
        break
      seq = str(ref[j].seq[position-5 : position+6]) #extract 11mer with A in middle
      if recs[i]['strand'] == '-':
        seq = seq_tools.reverse_complement(seq)
      kmers.append(seq)
      break

#remove errant kmers that do not have A in middle position
fixed_kmers = [i for i in kmers if not (i[5] != 'A')]


#randomly select seqs from kmer list to use as negative controls, trim seq to desired length and m6A position

#if using a list of different positive target motifs to avoid in controls 
if args.target_list != None:
  targets = []
  for item in args.target_list:
    temp_list = seq_tools.ambig_seq(item)
    for things in temp_list:
      targets.append(things)

#if using a single positive motif to avoid in controls
if args.target != None:   
  targets = seq_tools.ambig_seq(args.target) #specify a positive target motif to avoid 

#if not specifying any targets to avoid
if (args.target == None) and (args.target_list == None):
  targets = []

controls = []
while len(controls) < args.n:
  seq = random.choice(fixed_kmers)
  adjusted_seq = seq[6-(args.mod_position):6+(args.len - args.mod_position)] #trim
  if (adjusted_seq in targets) or (adjusted_seq in controls): #to avoid duplicate or match to positive target
    continue
  else:
    controls.append(adjusted_seq)
for i in controls:
  print(i)

