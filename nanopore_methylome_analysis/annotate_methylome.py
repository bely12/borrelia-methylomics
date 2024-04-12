##### required packages #####
import Bio
from Bio import SeqIO
from Bio import SeqRecord
from Bio import GenBank
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter



##### ARGUMENTS #####
parser = argparse.ArgumentParser(
    description='tallies up modified base counts within all coding and non-coding regions of the genome'
    ' \n', formatter_class=RawTextHelpFormatter)
parser.add_argument('-gbff', help='gbff annotation file, needs to match reference genome reads were mapped to\n')
parser.add_argument('-sites', help='tsv file containing modified base sites; use output from mod_mapper.py for proper formatting')
parser.add_argument('-out', help='prefix for output file')
args = parser.parse_args()


##### EXECUTION #####
### create table of cds regions from gbff file ###
gb_file = args.gbff
output = []
for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank"):
  for i in range(len(gb_record.features)):
    if gb_record.features[i].type == 'CDS':
      output.append({'element_id': gb_record.id, 'element_len': len(gb_record.seq), 'gene_id': gb_record.features[i].qualifiers['locus_tag'], 'start': gb_record.features[i].location.start, 'end': gb_record.features[i].location.end})
df = pd.DataFrame(output)


### add intergenic regions to table ###
regions = []
for plasmid in df.element_id.unique():
  df2 = df[(df['element_id'] == plasmid)].reset_index()
  regions.append({'seqid': plasmid, 'seq_len': df2.loc[0]['element_len'], 'start': 1, 'end': df2.loc[0]['start']-1, 'gene_id': 'intergenic'})
  for row in range(len(df2.index)-1):
    regions.append({'seqid': plasmid, 'seq_len': df2.loc[0]['element_len'], 'start': df2.loc[row]['start'], 'end': df2.loc[row]['end'], 'gene_id': df2.loc[row]['gene_id']})
    if row < len(df2.index):
      if df2.loc[row]['end']+1 < df2.loc[row+1]['start']:
        regions.append({'seqid': plasmid, 'seq_len': df2.loc[0]['element_len'], 'start': df2.loc[row]['end']+1, 'end': df2.loc[row+1]['start']-1, 'gene_id': 'intergenic'})
  regions.append({'seqid': plasmid, 'seq_len': df2.loc[0]['element_len'], 'start': df2.loc[row+1]['start'], 'end': df2.loc[row+1]['end'], 'gene_id': df2.loc[row+1]['gene_id']})
  regions.append({'seqid': plasmid, 'seq_len': df2.loc[0]['element_len'], 'start': df2.iloc[-1]['end']+1, 'end': df2.loc[row]['element_len'], 'gene_id': 'intergenic'})


### import and format mod sites table ###
mods = pd.read_table(args.sites, header=None)
mods = mods.rename(columns={0: 'motif', 1: 'plasmid_id', 2: 'position', 3: 'strand'})
for column in mods:
  if column != 'position':
    mods[column] = mods[column].str.strip() #remove leading and trailing spaces, not sure why this happens
mods = mods.to_dict('records')


### annotate the modified bases ###
results = []
for i in range(len(regions)):
  id = regions[i]['seqid']
  x = list(filter(lambda d: d['plasmid_id'] in id, mods))
  n = 0
  for j in range(len(x)):
    if ( x[j]['position'] >= regions[i]['start'] ) and ( x[j]['position'] <= regions[i]['end'] ):
      n+=1
  results.append(dict(regions[i], counts = n))


##### OUTPUT #####
# convert results into pandas df
annotated_mods = pd.DataFrame.from_dict(results)
# save to tsv file
annotated_mods.to_csv(args.out+'_methylome_annotations.tsv', sep='\t', index=False, header=True)