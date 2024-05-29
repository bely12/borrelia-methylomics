#!/usr/bin/env python

import Bio
from Bio import SeqIO
import re
import pandas as pd
import numpy as np
import random
from scipy import stats
from scipy.stats import ttest_1samp
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from tools_and_utilities import seq_tools
import argparse
from argparse import RawTextHelpFormatter

########## ARGUMENTS ##########
parser = argparse.ArgumentParser(
    description='This package contains tools for testing motif enrichment and avoidance hypotheses.\n'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-fasta', help='fasta file of cds\n'
                    ' \n')

parser.add_argument('-sim', help='fasta file of simulated cds\n'
                    ' \n')

parser.add_argument('-motif', help='file containing 1 motif per line.\n'
                     ' \n')

parser.add_argument('-tid', '--ncbi_tax_id', help='NCBI taxon id. Required. Default = 139 (Borrelia burgdorferi)\n'
                    ' \n', default = 139, type = int)

parser.add_argument('-out', help='pre-fix for output file', default=None)

args = parser.parse_args()


#define the motifs to search for
with open(args.motif, 'r') as motifs:
  motif_list = motifs.readlines()
for i in range(len(motif_list)):
  motif_list[i] = motif_list[i].replace("\n", "")


##### count motifs in real seq #####

#count motif occurences
individual_counts = seq_tools.motif_counter(motif_list, args.fasta)

#define multi fasta
seqList = list(SeqIO.parse(args.fasta, "fasta"))

#get a list of ids from multi fasta
id_list = []
for seq in seqList:
  id_list.append(seq.id)

#aggregate individual counts to totals for each consensus motif
final_list = []
for motif in motif_list:
  filtered_consensus = [x for x in individual_counts if x['consensus'] == motif]

  for gene in id_list:
    filtered_gene = [x for x in filtered_consensus if x['id'] == gene]

    final_list.append({'motif': motif,
                       'id': gene,
                       'plus': sum(item.get('fwd_match') for item in filtered_gene),
                       'minus': sum(item.get('rev_match') for item in filtered_gene),
                       'total_count': sum(item.get('total_match') for item in filtered_gene)})

df = pd.DataFrame(final_list)


##### count motifs in simulated seqs #####

#count motif occurences
sim_individual_counts = seq_tools.motif_counter(motif_list, args.sim)

#define multi fasta
seqList = list(SeqIO.parse(args.sim, "fasta"))

#get a list of ids from multi fasta
id_list = []
for seq in seqList:
  id_list.append(seq.id)

#aggregate individual counts to totals for each consensus motif
sim_final_list = []
for motif in motif_list:
  filtered_consensus = [x for x in sim_individual_counts if x['consensus'] == motif]

  for gene in id_list:
    filtered_gene = [x for x in filtered_consensus if x['id'] == gene]

    sim_final_list.append({'motif': motif,
                           'id': gene,
                           'plus': sum(item.get('fwd_match') for item in filtered_gene),
                           'minus': sum(item.get('rev_match') for item in filtered_gene),
                           'total_count': sum(item.get('total_match') for item in filtered_gene)})

df2 = pd.DataFrame(sim_final_list)
df2['id'] = df2['id'].str.replace(r'_sim.', '', regex=True)


##### 1 sample t test #####

#t-test
tstats = []
for seq_id in df.id.unique():
  #print(seq_id)
  temp_obs = df[df['id'] == seq_id]
  temp_sim = df2[df2['id'] == seq_id]
  for motif in df.motif.unique():
    #print(motif)
    temp_obs2 = temp_obs[temp_obs['motif'] == motif]
    temp_sim2 = temp_sim[temp_sim['motif'] == motif]
    #print(temp_sim2)
    #print('next')
    x = temp_obs2['total_count'].mean()
    y = temp_sim2['total_count'].mean()
    a = stats.ttest_1samp(temp_sim2['total_count'],x)
    if x == 0:
      x = 0.1
    if y == 0:
      y = 0.1
    #tstats.append({'id': seq_id,'motif': motif, 'test_stat': round(a[0],2), 'p_val': round(a[1],5), 'observed': x, 'sim_mean': round(y,2),'log2fc': np.round(np.log2(x/y),3)})
    tstats.append({'id': seq_id,'motif': motif, 'test_stat': round(a[0],2), 'p_val': float("{:.2e}".format(a[1])), 'observed': x, 'sim_mean': round(y,2),'log2fc': np.round(np.log2(x/y),3)})

for item in tstats:
  #if item['p_val'] < 0.005:
  print(item['id'],'\t',item['motif'],"\t",item['test_stat'],"\t",item['p_val'],'\t', item['observed'],'\t', item['sim_mean'],'\t', item['log2fc'])

if args.out != None:
    df3 = pd.DataFrame(tstats)
    df3['neg_log10'] = round(-np.log10(df3['p_val']),5)
    df3.to_csv(args.out+'_t-test_results.tsv', sep='\t', index=False, header=True)
