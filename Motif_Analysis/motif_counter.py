#!/usr/bin/env python

import Bio
from Bio import SeqIO
import seq_tools
import argparse
from argparse import RawTextHelpFormatter

#ARGUMENTS
parser = argparse.ArgumentParser(
    description='This package contains tools for testing motif enrichment and avoidance hypotheses.\n'
    ' \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('-fasta', help='fasta file of cds\n')
parser.add_argument('-motif', help='file containing 1 motif per line.\n')
parser.add_argument('-out', help='pre-fix for output file', default=None)
args = parser.parse_args()



#define the motifs to search for
with open(args.motif, 'r') as motifs:
  motif_list = motifs.readlines()
for i in range(len(motif_list)):
  motif_list[i] = motif_list[i].replace("\n", "")

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

#df = pd.DataFrame(final_list)
for item in final_list:
  print(item['motif'],'\t',
        item['id'],'\t',
        item['plus'],'\t',
        item['minus'],'\t',
        item['total_count'])

