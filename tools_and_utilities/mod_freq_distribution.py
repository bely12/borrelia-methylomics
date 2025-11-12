########## packages 
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from argparse import RawTextHelpFormatter

########## ARGUMENTS
parser = argparse.ArgumentParser(
    description='distribution of mod frequency values for all eligible positions'
    ' \n', formatter_class=RawTextHelpFormatter)


parser.add_argument('-bed', help='tsv bed file with mod calls and kmers. kmers must be in the 3rd column of the table, no header\n')
parser.add_argument('-min_cov', type=int, default=10, help='mininum position coverage to filter for')
parser.add_argument('-out', default=None, help='prefix for output file')
args = parser.parse_args()


#prepare bed file and convert to list of row dictionaries 
bed = pd.read_table(args.bed, sep = '\t', header=None)

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

# filter for specified coverage
filtered_bed = bed[bed['cov'] >= args.min_cov]

if filtered_bed.empty:
    raise ValueError(f"No positions with coverage >= {args.min_cov}")

# histogram
plt.figure(figsize=(7,4))
plt.hist(filtered_bed['mod_freq'], bins=100, color='blue', edgecolor='black')
plt.xlabel('mod frequency')
plt.ylabel('count')

out_plot = f"{args.out}.png"
plt.tight_layout()
plt.savefig(out_plot, dpi=300)