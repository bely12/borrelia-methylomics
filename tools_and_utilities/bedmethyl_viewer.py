# its simple - look at, sort, get info from bedmethyl file
########## packages 
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter

##### Arguments #####

parser = argparse.ArgumentParser(
    description='View and filter a bedmethyl file produced by the Modkit [pileup] output'
    ' \n', formatter_class=RawTextHelpFormatter)

### Upload inputs (bedmethyl file from modkit)
parser.add_argument('-bed', help='tsv bed file with mod calls and kmers. kmers must be in the 3rd column of the table, no header\n')

### filtering options
parser.add_argument('-depth', default=None, type=int, help='minimum sequencing depth to use, inclusive')
parser.add_argument('-mod_freq', default=None, type=float, help='percentage of reads that need to be predicted as modified for that position to call position 6mA')

### strand summary
parser.add_argument('-strand', const=True, nargs='?')

### viewing and sorting
parser.add_argument('-view', const=True, nargs='?')
parser.add_argument('-n',type=int,default=10, help='number of rows to display')
parser.add_argument('-sort',default=None, help='what column to sort by')

args = parser.parse_args()


##### import bed #####
bed = pd.read_table(args.bed, sep = '\t', header=None)
bed[['cov', 'percent_mod', 'N_mod', 'N_can', 'N_other', 'N_del', 'N_fail', 'N_diff', 'N_nocall']] = bed[9].str.split(' ', expand=True)
bed = bed.rename(columns={0: 'chromosome', 1: 'start', 2: 'end', 3: 'mod', 4: 'score', 5: 'strand'})
bed = bed.drop(columns=[6,7,8,9])
bed['percent_mod'] = bed['percent_mod'].astype(float)
bed['cov'] = bed['cov'].astype(int)

##### filtering #####
if args.depth != None:
  bed = bed[ bed['cov'] >= args.depth]

if args.mod_freq != None:
  bed = bed[ bed['percent_mod'] >= args.mod_freq]

##### sorting #####
if args.sort != None:
  bed = bed.sort_values(by=[args.sort], ascending=False)

##### viewing #####
if args.view == True:
  table = bed.iloc[:args.n] 
  print(table)
  print('Total # of positions:', len(bed))

##### strand info #####
if args.strand == True:
  counts = bed['strand'].value_counts()
  print('positions per strand:')
  for i in range(len(counts.index)):
    print(counts.index[i], counts.values[i])

