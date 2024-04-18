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

  #get names of plasmids --- this may not work for all gbff files (and definitely not for non borrelia genomes). Add {'name': name} to output if using this
  # if 'plasmid' in gb_record.description:
  #   name = ''.join(re.findall(r'lp.*,|cp.*,', gb_record.description))
  #   name = name.replace(',','')
  # else:
  #   name = 'chromosome'

  #parse features
  for i in range(len(gb_record.features)):

    #get all cds records
    if gb_record.features[i].type == 'CDS':

      #get gene names, try for old locus, if not there, get new locus
      if gb_record.features[i].qualifiers.get('old_locus_tag') is not None:
        hold_gene_id = ''.join(gb_record.features[i].qualifiers['old_locus_tag'])
      else:
        hold_gene_id = ''.join(gb_record.features[i].qualifiers['locus_tag'])

      #get protein ids, if they are there. If not, assign na
      if gb_record.features[i].qualifiers.get('protein_id') is not None:
        hold_prot_id = ''.join(gb_record.features[i].qualifiers['protein_id'])
      else:
        hold_prot_id = 'na'

      #get strand orientation for cds
      if gb_record.features[i].location.strand == 1:
        hold_strand = '+'
      else:
        hold_strand = '-'

      #store the needed info as a list of dictionaries
      output.append({'element_id': gb_record.id,
                    'element_len': len(gb_record.seq),
                    'gene_id': hold_gene_id,
                    'protein_id': hold_prot_id,
                    'start': gb_record.features[i].location.start,
                    'end': gb_record.features[i].location.end,
                    'strand': hold_strand})
df = pd.DataFrame(output)


### add intergenic regions to table ###
regions = []
for plasmid in df.element_id.unique():
  df2 = df[(df['element_id'] == plasmid)].reset_index()
  regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': 1, 'end': df2.loc[0]['start']-1, 'gene_id': 'intergenic', 'protein_id': 'na', 'strand': 'na'})
  for row in range(len(df2.index)-1):
    regions.append({'element_id': plasmid,'len': df2.loc[0]['element_len'], 'start': df2.loc[row]['start'], 'end': df2.loc[row]['end'], 'gene_id': df2.loc[row]['gene_id'], 'protein_id': df2.loc[row]['protein_id'], 'strand': df2.loc[row]['strand']})
    if row < len(df2.index):
      if df2.loc[row]['end']+1 < df2.loc[row+1]['start']:
        regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': df2.loc[row]['end']+1, 'end': df2.loc[row+1]['start']-1, 'gene_id': 'intergenic', 'protein_id': 'na', 'strand': 'na'})
  regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': df2.loc[row+1]['start'], 'end': df2.loc[row+1]['end'], 'gene_id': df2.loc[row+1]['gene_id'], 'protein_id': df2.loc[row+1]['protein_id'], 'strand': df2.loc[row+1]['strand']})
  regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': df2.iloc[-1]['end']+1, 'end': df2.loc[row]['element_len'], 'gene_id': 'intergenic','protein_id': 'na', 'strand': 'na'})


### import and format mod sites table ###
mods = pd.read_table(args.sites, header=None)
mods = mods.rename(columns={0: 'motif', 1: 'plasmid_id', 2: 'position', 3: 'strand', 4: 'status'})
for column in mods:
  if column != 'position':
    mods[column] = mods[column].str.strip() #remove leading and trailing spaces, not sure why the spaces are there
mods = mods.to_dict('records')


### annotate the modified bases ###
results = []
for i in range(len(regions)):
  
  #filter to search faster
  id = regions[i]['element_id']
  x = list(filter(lambda d: d['plasmid_id'] in id, mods))
  
  #counters for mods (target mod, non-mod, and off target mods)
  mod = 0
  non_mod = 0
  off_target_mod = 0 
  
  #search and count
  for j in range(len(x)):
    if ( x[j]['position'] >= regions[i]['start'] ) and ( x[j]['position'] <= regions[i]['end'] ):
      if x[j]['status'] == 'modified' and x[j]['motif'] != 'off_target':
        mod+=1
      if x[j]['status'] == 'unmodified' and x[j]['motif'] != 'off_target':
        non_mod+=1
      if x[j]['motif'] == 'off_target':
        off_target_mod+=1
  
  #add counts to table 
  results.append(dict(regions[i], target_mod = mod, target_no_mod = non_mod, off_target = off_target_mod))


##### OUTPUT #####

# convert results into pandas df
annotated_mods = pd.DataFrame.from_dict(results)
annotated_mods['feature_len'] = annotated_mods['end'] - annotated_mods['start']
annotated_mods['mod_rate'] = annotated_mods['target_mod'] / annotated_mods['feature_len']
annotated_mods = annotated_mods[['element_id', 'len', 'gene_id', 'protein_id', 'start', 'end', 'feature_len', 'strand', 'target_mod', 'target_no_mod', 'off_target','mod_rate']]
# save to tsv file
annotated_mods.to_csv(args.out+'_methylome_annotations.tsv', sep='\t', index=False, header=True)