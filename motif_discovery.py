########## packages 
import Bio
from Bio import SeqIO
import pandas as pd
import numpy as np
import random
import seq_tools
import argparse
from argparse import RawTextHelpFormatter

########## ARGUMENTS
parser = argparse.ArgumentParser(
    description=''
    ' \n', formatter_class=RawTextHelpFormatter)


##### arguments for making kmer fasta files for motif enrichment analysis
parser.add_argument('-kmer_fasta', const=True, nargs='?')
parser.add_argument('-bed', help='tsv bed file with mod calls and kmers. kmers must be in the 3rd column of the table, no header\n')
parser.add_argument('-ref', default = None, help='reference genome in fasta or multi fasta format, first seq will be used to generate control seqs')
parser.add_argument('-nt_convert', default= True, help='converts M to A in kmers')
parser.add_argument('-out',default=None, help='prefix for output fasta file names')

##### arguments for reading position weight matrix and generating consensus motif and all possible sequences
parser.add_argument('-streme_out', help='probability matrix txt file from streme output')
parser.add_argument('-motif_num', type=int, help='which  number motif from streme_out you want to use')
parser.add_argument('-get_consensus', default= False, help='takes position weight matrix and returns consensus motif wtih all possible seqs' )
parser.add_argument('-threshold',type=float, default=0.85, help='freq threshold for calling a symbol')

##### arguments for calculating the total percent modification across the genome for a motif of interest
parser.add_argument('-percent_mod', const=True, nargs='?')
parser.add_argument('-motif', default=None, help='motif of interst')
parser.add_argument('-mod_pos', type=int, help='position of predicted mod base in motif')
parser.add_argument('-mod_call_thresh', type=float, default=0.5, help='threshold for calling a position modified in bed file')
parser.add_argument('-min_cov', type=int, default=1, help='mininum position coverage to filter for')

args = parser.parse_args()

############# get positive and control fasta files for motif enrichment analysis 

if args.kmer_fasta == True:
  #upload bed file from mCaller
  bed = pd.read_table(args.bed, sep = '\t', header=None)
  #bed

  #prepare 11mer sequences
  if args.nt_convert == True:
    seqs = bed[3].tolist()
    #seqs
    mod_seqs = []
    for i in range(len(seqs)):
      mod_seqs.append(seqs[i].replace('M','A'))
    #mod_seqs
  else:
    mod_seqs = bed[3].tolist()

  #create multi fasta file for 11-mers surrounding methylated position
  k = 1
  for seq in mod_seqs:
      print('>motif_',k,'\n',seq, sep='', file = open(args.out+'_mod_kmers.fasta', 'a'))
      k += 1

  #create multi fasta for control seeks
  if args.ref != None:
    ref_fasta = list(SeqIO.parse(args.ref, "fasta"))
    window = len(mod_seqs[0]) #length of control seek
    k = 1
    for i in range(len(mod_seqs)): 
      position = random.randint(200, len(ref_fasta[0].seq)-200)
      seq = ref_fasta[0].seq[position:position+window]
      print('>control_', k,'\n', seq, sep='', file = open(args.out+'_control_kmers.fasta', 'a'))
      k += 1


############# Module: input the streme.txt output, return a pos weight matrix specified motif or get a consensus seq

if args.streme_out != None:
  matrix = seq_tools.get_weight_matrix(args.streme_out, args.motif_num)
  matrix.index = matrix.index + 1
  
  if args.get_consensus == False:
    print(matrix)
  
  else:
    consensus_motif = seq_tools.get_consensus_motif(matrix, args.threshold)
    print('Consenus motif: ',consensus_motif['consensus'])
    #print('\n')
    print('Possible seqs:')
    for i in consensus_motif['seqs']:
      print(i)


############# input an amibigous seq, get back a list of all possible seqs

if args.motif != None and args.percent_mod != True:
  print('Input seq: ',args.motif)
  print('Possible seqs:')
  all_seqs = seq_tools.ambig_seq(args.motif)
  for i in range(len(all_seqs)):
    print(all_seqs[i])


############# calculate total percent methylated for specified motif

if args.percent_mod == True:
  #prepare bed file
  bed = pd.read_table(args.bed, sep = '\t', header=None)
  bed = bed.rename(columns={0: 'chromosome', 1: 'start', 2: 'end', 3: 'old_motif', 4: 'mod_freq', 5: 'strand', 6: 'cov'})
  #drop bad columns (drops a column if it doesn't have a kmer, appears to just be the last one)
  for i in range(len(bed.index)):
    if type(bed['old_motif'][i]) != str:
      bed.drop(index = i, inplace=True)
  #len(bed.index)
  if args.nt_convert == True:
    seqs = bed['old_motif'].tolist()
    #seqs
    mod_seqs = []
    for i in range(len(seqs)):
      mod_seqs.append(seqs[i].replace('M','A'))
    #mod_seqs
  else:
    mod_seqs = bed[3].tolist()
  #add fixed seqs as a column in table
  bed['kmer'] = mod_seqs
  bed = bed.drop(['old_motif'], axis=1)
  #convert bed df to dictionary
  bed2 = bed.to_dict('records')
  
  #generate list of motif possibilities, find all occurances, and calculate percent methylated for each
  motif_list = seq_tools.ambig_seq(args.motif)
  results = seq_tools.get_mod_frequncy(motif_list, bed2, mod_pos = args.mod_pos, mod_call_thresh = args.mod_call_thresh)
  #summarize results
  print('motif','\t','occurences','\t','mod','\t','percent')
  for i in range(len(results)):
    print(results[i]['motif'],'\t',results[i]['occurences'],'\t',results[i]['modified'],'\t',results[i]['percent_mod'])
  print('Total percent modfied for all motifs = ',
        round(sum(item.get('modified') for item in results) / sum(item.get('occurences') for item in results) * 100, 2))
  
  # all_occur = []
  # all_mod = []
  # for i in range(len(results)):
  #   all_occur.append(results[i]['occurences'])
  #   all_mod.append(results[i]['modified'])
  # print('Total percent modfied for all motifs = ',round((np.sum(all_mod) / np.sum(all_occur))*100,2))
    
