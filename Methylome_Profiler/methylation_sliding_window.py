
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter

#description
parser = argparse.ArgumentParser(
    description='This program does a sliding window analysis for methylated sites in a genome.\n'
    ' \n', formatter_class=RawTextHelpFormatter)

#arguments
parser.add_argument('-sites', help='file containing methylation sites; use output from mod_mapper.py\n')
parser.add_argument('-seq_len', help='file containing a column with seq id and a column with seq length\n')
parser.add_argument('-win', type=int, help='size of window\n')
parser.add_argument('-step', type=int, help='size of step for slide\n')
parser.add_argument('-all', const=True, nargs='?', help='returns both methylated target sites and off-target sites')
parser.add_argument('-out', help='pre-fix for output file', default=None)
args = parser.parse_args()


#sliding window function 
def sliding_window(mod_dict, win_size, step):
  #initialize empty lists to collect info
  output = []
  id_set = []
  #get a list of unique sequence id's to filter with
  for rec in mod_dict:
    if rec['id'] not in id_set:
      id_set.append(rec['id'])
  #filter for a sequence id
  for id in id_set:
    temp_dict = list(filter(lambda x: x['id'] == id, mod_dict))
    #count number of mods in window
    for i in range(0, temp_dict[0]['len'], step):
      #break loop if out of range
      if i > temp_dict[0]['len'] - win_size:
        break
      #filter by position and use len of filtered list as the mod count
      else:
        mod_count = len(list(filter(lambda x: (x['pos'] >= i) & (x['pos'] < i + win_size), temp_dict)))
        output.append({'id': id, 'start': i, 'end': i + win_size, 'mods': mod_count})
  #return a new list of dictionaries with mod counts for each window
  return output


#load sites table and name columns
df = pd.read_table(args.sites, header=None, sep='\t')
df = df.rename(columns={0: 'motif', 1: 'id', 2: 'pos', 3: 'strand', 4: 'status'})

#strip leading/lagging spaces from strings
for column in df:
  if column != 'pos':
    df[column] = df[column].str.strip()

#filter for types of methylation (target site only, or off-target site as well)
if args.all != None:
  df = df.loc[ df['status'] == 'modified' ]
else:
  df = df.loc[(df['motif'] != 'off_target') & (df['status'] == 'modified')]

mod_dict = df.to_dict('records')

#load table with seq id's and lengths 
df2 = pd.read_table(args.seq_len, sep='\t')
seq_len = df2.to_dict('records')

#add seq lengths to dictionaries 
out = []
for i in range(len(seq_len)):
  for j in range(len(mod_dict)):
    if seq_len[i]['id'] == mod_dict[j]['id']:
      out.append(dict(mod_dict[j], len = seq_len[i]['len']))
mod_dict = out


#sliding window analysis 
results = sliding_window(mod_dict, args.win, args.step)

#save results to file
if args.out != None:
  df3 = pd.DataFrame(results)
  df3.to_csv(args.out+'_windows.tsv', sep='\t', index=False, header=True)

#print results 
else:
  for rec in results: 
    print(rec['id'],'\t',rec['start'],'\t',rec['end'],'\t',rec['mods'])
    

