########## packages 
import Bio
from Bio import SeqIO
import pandas as pd
import numpy as np
import seq_tools
from scipy.stats import binomtest
import argparse
from argparse import RawTextHelpFormatter

########## ARGUMENTS
parser = argparse.ArgumentParser(
    description='calculates what percentage of a specified motif is modified in the genome'
    ' \n', formatter_class=RawTextHelpFormatter)

#required
parser.add_argument('-bed', help='tsv bed file with mod calls and kmers. kmers must be in the 3rd column of the table, no header\n')
parser.add_argument('-ref', default = None, help='reference genome in fasta or multi fasta format, first seq will be used to generate control seqs')
parser.add_argument('-motif', default=None, help='motif of interst')
parser.add_argument('-control', default=None, help='list of seqs from txt file, one seq per line')
parser.add_argument('-mod_pos', type=int, help='position of predicted mod base in motif')

#args with defaults that can be customized
parser.add_argument('-mod_call_thresh', type=float, default=50.0, help='threshold for calling a position modified in bed file')
parser.add_argument('-min_cov', type=int, default=10, help='mininum position coverage to filter for')

#save percent modfiied table with binomial distribution test as output tsv file
parser.add_argument('-out', default=None, help='prefix for output file name')

args = parser.parse_args()


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


# define negative controls
with open(args.control, 'r') as motifs:
  motif_list = motifs.readlines()
for i in range(len(motif_list)):
  motif_list[i] = motif_list[i].replace("\n", "")

# get percent modified for each motif in negative control list and calculate mean for binomial distribution test
results = seq_tools.get_mod_frequency(motif_list, bed3, mod_pos = args.mod_pos, mod_call_thresh = args.mod_call_thresh, min_cov = args.min_cov)
total_occurences = sum(item.get('occurences') for item in results)
total_modified = sum(item.get('modified') for item in results)
prob = round((total_modified / total_occurences) * 100, 2)


### pecent modified and binomial distribution test on each candidate sequence

# define candidate motif
motif_list = seq_tools.ambig_seq(args.motif)
results = seq_tools.get_mod_frequncy(motif_list, bed3, mod_pos = args.mod_pos, mod_call_thresh = args.mod_call_thresh, min_cov = args.min_cov)

# make column names for results
print('motif','\t','occurences','\t','mod','\t','percent','\t','pvalue','\t','95% CI')

for i in range(len(results)):
  binom_test = binomtest(k=results[i]['modified'], n=results[i]['occurences'], p=prob, alternative='greater')
  ci = binom_test.proportion_ci(confidence_level = 0.95)
  standev = round((ci[1]-ci[0])/2,2)
  ci_range = round(ci[0] + standev, 2)

  print(results[i]['motif'],'\t',
        results[i]['occurences'],'\t',
        results[i]['modified'],'\t',
        results[i]['percent_mod'], '\t',
        round(binom_test.pvalue,4),'\t',
        ci_range,'+/-',standev)
  
  # save results as an output file
  if args.out!=None:
    print(results[i]['motif'],'\t',
          results[i]['occurences'],'\t',
          results[i]['modified'],'\t',
          results[i]['percent_mod'], '\t',
          round(binom_test.pvalue,4),'\t',
          ci_range,'+/-',standev,
          file=open(args.out+'_mod_percentage_results.tsv', 'a'))

### percent modified and binomial distribution test for consensus motif
total_occurences = sum(item.get('occurences') for item in results)
total_modified = sum(item.get('modified') for item in results)
binom_test = binomtest(k=total_modified, n=total_occurences, p=prob, alternative='greater')
ci = binom_test.proportion_ci(confidence_level = 0.95)
standev = round((ci[1]-ci[0])/2,2)
ci_range = round(ci[0] + standev, 2)

#print results
print('\n')
print('Total percent modfied for all motifs = ',round((total_modified / total_occurences) * 100, 2))
print('Binomial sampling test:','\n',
      'k=',total_modified,'\n',
      'n=',total_occurences,'\n',
      'prob=',prob,'\n',
      'p value =',round(binom_test.pvalue,4),'\n', '95% CI =', ci_range,'+/-',standev)
