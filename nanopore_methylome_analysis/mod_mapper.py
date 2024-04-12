##### packages #####
import Bio
from Bio import SeqIO
import pandas as pd
import seq_tools
import argparse
from argparse import RawTextHelpFormatter

##### ARGUMENTS #####
parser = argparse.ArgumentParser(
    description='gives coordinates of all likely modified adenines within identified RM target sites'
    ' \n', formatter_class=RawTextHelpFormatter)


parser.add_argument('-bed', help='tsv bed file with mod calls, no header\n')
parser.add_argument('-ref', default = None, help='reference genome in fasta or multi fasta format, first seq will be used to generate control seqs')
parser.add_argument('-motif', default=None, help='motif of interst')
parser.add_argument('-motif_file', default=None, help='list of seqs from txt file, one seq per line')
parser.add_argument('-mod_pos', type=int, help='position of predicted mod base in motif')
parser.add_argument('-mod_call_thresh', type=float, default=50.0, help='threshold for calling a position modified in bed file')
parser.add_argument('-min_cov', type=int, default=10, help='mininum position coverage to filter for')
args = parser.parse_args()

##### EXECUTION #####
#import reference genome and change all letters to uppercase
ref = list(SeqIO.parse(args.ref, "fasta"))
for i in range(len(ref)):
  ref[i].seq = ref[i].seq.upper()

#prepare bed file and convert to list of row dictionaries 
bed = pd.read_table(args.bed, sep = '\t', header=None)
bed[['cov', 'mod_freq', 'N_mod', 'N_can', 'N_other', 'N_del', 'N_fail', 'N_diff', 'N_nocall']] = bed[9].str.split(' ', expand=True)
bed = bed.rename(columns={0: 'chromosome', 1: 'start', 2: 'end', 3: 'mod', 4: 'score', 5: 'strand'})
bed = bed.drop(columns=[6,7,8,9])
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
      seq = str(ref[j].seq[position-5 : position+6])
      if recs[i]['strand'] == '-':
        seq = seq_tools.reverse_complement(seq)
      #make a new list of dictionaries by adding 11mer as key value pair to existing bedmethyl dictionaries
      bed2.append(dict(recs[i], kmer = seq))
      break
#bed2

#remove errant kmers (don't have potential 6mA in correct position)
bed3 = [i for i in bed2 if not ( (len(i['kmer']) != 11) or (i['kmer'][5] != 'A') )]

#if input is a single motif with the -motif arg 
if args.motif != None:
  motif_list = seq_tools.ambig_seq(args.motif)

#if input is a txt file containing a list of motifs, 1 per line (no ambiguous NT's) *use with neg control seqs
if args.motif_file != None:
  with open(args.motif_file, 'r') as motifs:
    motif_list = motifs.readlines()
  for i in range(len(motif_list)):
    motif_list[i] = motif_list[i].replace("\n", "")

#get coordinates of modified adenines 
results = seq_tools.mod_mapper(motif_list, bed3, mod_pos = args.mod_pos, mod_call_thresh = args.mod_call_thresh, min_cov = args.min_cov)

##### OUTPUT #####
#print results
for i in range(len(results)):
  print(results[i]['motif'],'\t',results[i]['plasmid_id'],'\t',results[i]['position'],'\t',results[i]['strand'])