########## packages 
import Bio
from Bio import SeqIO
import pandas as pd
import numpy as np
from scipy.stats import binomtest
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from tools_and_utilities import seq_tools
from tools_and_utilities import mod_tools
import argparse
from argparse import RawTextHelpFormatter

########## ARGUMENTS
parser = argparse.ArgumentParser(
    description='calculates what percentage of a specified motif is modified in the genome'
    ' \n', formatter_class=RawTextHelpFormatter)


parser.add_argument('-bed', help='tsv bed file with mod calls and kmers. kmers must be in the 3rd column of the table, no header\n')
parser.add_argument('-ref', default = None, help='reference genome in fasta or multi fasta format, first seq will be used to generate control seqs')
parser.add_argument('-motif', default=None, help='motif of interst')
parser.add_argument('-motif_file', default=None, help='list of seqs from txt file, one seq per line')
parser.add_argument('-mod_base', default = 'A', help='what modified base ex. A or C')
parser.add_argument('-mod_pos', type=int, help='position of predicted mod base in motif')
parser.add_argument('-mod_call_thresh', type=float, default=50.0, help='threshold for calling a position modified in bed file')
parser.add_argument('-min_cov', type=int, default=10, help='mininum position coverage to filter for')
parser.add_argument('-out', default=None, help='prefix for output file')
args = parser.parse_args()


#import reference genome and change all letters to uppercase
ref = list(SeqIO.parse(args.ref, "fasta"))
for i in range(len(ref)):
  ref[i].seq = ref[i].seq.upper()

#prepare bed file and convert to list of row dictionaries 
bed = pd.read_table(args.bed, sep = '\t', header=None)
#bed[['cov', 'mod_freq', 'N_mod', 'N_can', 'N_other', 'N_del', 'N_fail', 'N_diff', 'N_nocall']] = bed[9].str.split(' ', expand=True) #new modkit output separates columns better, so don't need this
#bed = bed.rename(columns={0: 'chromosome', 1: 'start', 2: 'end', 3: 'mod', 4: 'score', 5: 'strand'})
#bed = bed.drop(columns=[6,7,8,9]) #need this with old modkit output, use below for new output
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
      #seq = str(ref[j].seq[position-5 : position+6])
      seq = str(ref[j].seq[position-15 : position+16]) # try this for 31mers instead of 11mers for longer motifs
      if recs[i]['strand'] == '-':
        seq = seq_tools.reverse_complement(seq)
      #make a new list of dictionaries by adding 11mer as key value pair to existing bedmethyl dictionaries
      bed2.append(dict(recs[i], kmer = seq))
      break
#bed2

#remove errant kmers (don't have potential 6mA in correct position)
#bed3 = [i for i in bed2 if not ( (len(i['kmer']) != 11) or (i['kmer'][5] != args.mod_base) )]
bed3 = [i for i in bed2 if not ( (len(i['kmer']) != 31) or (i['kmer'][15] != args.mod_base) )]

#generate list of motif possibilities, find all occurances, and calculate percent methylated for each
#if input is a single motif with the --motif arg 
if args.motif != None:
  motif_list = seq_tools.ambig_seq(args.motif)

#if input is a txt file containing a list of motifs, 1 per line (no ambiguous NT's) *use with neg control seqs
if args.motif_file != None:
  with open(args.motif_file, 'r') as motifs:
    motif_list = motifs.readlines()
  for i in range(len(motif_list)):
    motif_list[i] = motif_list[i].replace("\n", "")

#get percent modified for each motif 
results = mod_tools.get_mod_frequency(motif_list, bed3, mod_pos = args.mod_pos, mod_call_thresh = args.mod_call_thresh, min_cov = args.min_cov)


### summarize individual results
#set column name for results table 
print('motif','\t','occurences','\t','mod','\t','percent')

#loop through results to perform statistical test and print data
for i in range(len(results)):
  print(results[i]['motif'],'\t',
            results[i]['occurences'],'\t',
            results[i]['modified'],'\t',
            results[i]['percent_mod'])

### summarize total results
total_occurences = sum(item.get('occurences') for item in results)
total_modified = sum(item.get('modified') for item in results)

#print total results info 
print('\n')
print('found in genome = ',total_occurences)
print('modified in genome = ', total_modified)
print('Total percent modfied for all motifs = ',
      round((total_modified / total_occurences) * 100, 2))

### output file ###
if args.out != None:
  df = pd.DataFrame.from_dict(results)
  df.to_csv(args.out+'_mod_calc.tsv', sep='\t', index=False, header=True)