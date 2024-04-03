########## packages 
import Bio
from Bio import SeqIO
import pandas as pd
import seq_tools
import argparse
from argparse import RawTextHelpFormatter

##### Arguments #####

parser = argparse.ArgumentParser(
    description=''
    ' \n', formatter_class=RawTextHelpFormatter)

### Choose one of these 3 (all kmers, likely modified kmers, or likely unmodified kmers)
parser.add_argument('-all', const=True, nargs='?')
parser.add_argument('-modified', const=True, nargs='?')
parser.add_argument('-unmodified', const=True, nargs='?')

### Upload inputs (bedmethyl file from modkit and reference genome)
parser.add_argument('-bed', help='tsv bed file with mod calls and kmers. kmers must be in the 3rd column of the table, no header\n')
parser.add_argument('-ref', default = None, help='reference genome in fasta or multi fasta format, first seq will be used to generate control seqs')

### filtering for minimum sequencing depth and criteria for calling a position 6mA
parser.add_argument('-depth', default=10, type=int, help='minimum sequencing depth to use, inclusive')
parser.add_argument('-mod_threshold', default=50.0, type=float, help='percentage of reads that need to be predicted as modified for that position to call position 6mA')

args = parser.parse_args()


##### Execution #####


### format input files

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



### sepcify what kmer set to extract

#get all kmers
if args.all == True:
  df = bed[ bed['cov'] >= args.depth] #filter for positions with at least 10X depth 
  recs = df.to_dict('records')

#get kmers for likely 6mA positions
if args.modified == True:
  df = bed[ (bed['percent_mod'] > args.mod_threshold) & (bed['cov'] >= args.depth)]
  recs = df.to_dict('records')

#get kmers for unlikely 6mA positions
if args.unmodified == True:
  df = bed[ (bed['percent_mod'] < args.mod_threshold) & (bed['cov'] >= args.depth)]
  recs = df.to_dict('records')



### extract 11mers for every adenine 

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



### print results in fasta format

k = 1
for seq in fixed_kmers:
  print('>motif_',k,'\n',seq, sep='')
  k += 1
