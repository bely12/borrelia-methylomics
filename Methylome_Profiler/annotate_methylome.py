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

# go through each plasmid or chromosome at a time
for plasmid in df.element_id.unique():
  df2 = df[(df['element_id'] == plasmid)].reset_index()
  
  # only add an intergenic row at the start if the first gene starts after position 1 in the sequence
  if df2.loc[0]['start'] > 1:
    regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': 1, 'end': df2.loc[0]['start']-1, 'gene_id': 'intergenic', 'protein_id': 'na', 'strand': 'na'})
  
  # if there is only one gene in the entire plasmid, add a row with that information, then move on to next plasmid
  if len(df2.index) == 1:
    regions.append({'element_id': plasmid,'len': df2.loc[0]['element_len'], 'start': df2.loc[0]['start'], 'end': df2.loc[0]['end'], 'gene_id': df2.loc[0]['gene_id'], 'protein_id': df2.loc[0]['protein_id'], 'strand': df2.loc[0]['strand']})

    # add a row for the intergenic space after the single gene, if it exists
    if df2.loc[0]['element_len'] > df2.loc[0]['end']:
      regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': df2.loc[0]['end']+1, 'end': df2.loc[0]['element_len'], 'gene_id': 'intergenic', 'protein_id': 'na', 'strand': 'na'})
  
  # if there is more than 1 gene in the plasmid, we will cycle through the rows one at a time but not get to the last one
  else:
    for row in range(len(df2.index)-1):

      # add row with gene info 
      regions.append({'element_id': plasmid,'len': df2.loc[0]['element_len'], 'start': df2.loc[row]['start'], 'end': df2.loc[row]['end'], 'gene_id': df2.loc[row]['gene_id'], 'protein_id': df2.loc[row]['protein_id'], 'strand': df2.loc[row]['strand']})
      
      # add a row for an intergenic region if it exists
      if df2.loc[row]['end']+1 < df2.loc[row+1]['start']:
        regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': df2.loc[row]['end']+1, 'end': df2.loc[row+1]['start']-1, 'gene_id': 'intergenic', 'protein_id': 'na', 'strand': 'na'})
    
    # add last row of gene info
    regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': df2.loc[row+1]['start'], 'end': df2.loc[row+1]['end'], 'gene_id': df2.loc[row+1]['gene_id'], 'protein_id': df2.loc[row+1]['protein_id'], 'strand': df2.loc[row+1]['strand']})
    
    # if the last gene on plasmid doesn't go to the end, add a final row for the intergenic space 
    if df2.loc[row+1]['element_len'] > df2.loc[row+1]['end']:
      regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': df2.iloc[-1]['end']+1, 'end': df2.loc[row]['element_len'], 'gene_id': 'intergenic','protein_id': 'na', 'strand': 'na'})

### this chunk makes a few minor errors that might mattter for some genomes and gbff files. Remove if the chunk above resolves problem!
# for plasmid in df.element_id.unique():
#   df2 = df[(df['element_id'] == plasmid)].reset_index()
#   regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': 1, 'end': df2.loc[0]['start']-1, 'gene_id': 'intergenic', 'protein_id': 'na', 'strand': 'na'})
#   for row in range(len(df2.index)-1):
#     regions.append({'element_id': plasmid,'len': df2.loc[0]['element_len'], 'start': df2.loc[row]['start'], 'end': df2.loc[row]['end'], 'gene_id': df2.loc[row]['gene_id'], 'protein_id': df2.loc[row]['protein_id'], 'strand': df2.loc[row]['strand']})
#     if row < len(df2.index):
#       if df2.loc[row]['end']+1 < df2.loc[row+1]['start']:
#         regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': df2.loc[row]['end']+1, 'end': df2.loc[row+1]['start']-1, 'gene_id': 'intergenic', 'protein_id': 'na', 'strand': 'na'})
#   regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': df2.loc[row+1]['start'], 'end': df2.loc[row+1]['end'], 'gene_id': df2.loc[row+1]['gene_id'], 'protein_id': df2.loc[row+1]['protein_id'], 'strand': df2.loc[row+1]['strand']})
#   regions.append({'element_id': plasmid, 'len': df2.loc[0]['element_len'], 'start': df2.iloc[-1]['end']+1, 'end': df2.loc[row]['element_len'], 'gene_id': 'intergenic','protein_id': 'na', 'strand': 'na'})


### import and format mod sites table ###
mods = pd.read_table(args.sites, header=None)
mods = mods.rename(columns={0: 'motif', 1: 'plasmid_id', 2: 'position', 3: 'strand', 4: 'status'})
for column in mods:
  if column != 'position':
    mods[column] = mods[column].str.strip() #remove leading and trailing spaces, not sure why the spaces are there
mods = mods.to_dict('records')


### annotate the modified bases ###
results = []
new_mod_sites = [] #*test

for i in range(len(regions)):
  
  #filter to search faster
  id = regions[i]['element_id']
  x = list(filter(lambda d: d['plasmid_id'] in id, mods))
  
  #counters
  mod = 0
  non_mod = 0
  off_target_mod = 0 
  
  #search and count
  for j in range(len(x)):
    if ( x[j]['position'] >= regions[i]['start'] ) and ( x[j]['position'] <= regions[i]['end'] ):
      if x[j]['status'] == 'modified' and x[j]['motif'] != 'off_target':
        mod+=1
        new_mod_sites.append(dict(x[j], annot = regions[i]['gene_id'])) #*test
      if x[j]['status'] == 'unmodified' and x[j]['motif'] != 'off_target':
        non_mod+=1
        new_mod_sites.append(dict(x[j], annot = regions[i]['gene_id'])) #*test
      if x[j]['motif'] == 'off_target':
        off_target_mod+=1
        new_mod_sites.append(dict(x[j], annot = regions[i]['gene_id'])) #*test
  
  #add counts to table 
  results.append(dict(regions[i], target_mod = mod, target_no_mod = non_mod, off_target = off_target_mod))


##### OUTPUT #####

# convert results into pandas df
annotated_mods = pd.DataFrame.from_dict(results)
annotated_mods['feature_len'] = (annotated_mods['end'] - annotated_mods['start']) +1
annotated_mods['mod_rate'] = round((annotated_mods['target_mod'] / annotated_mods['feature_len']) * 1000, 2)
annotated_mods['unmodified_rate'] = round((annotated_mods['target_no_mod'] / annotated_mods['feature_len']) * 1000, 2)
annotated_mods = annotated_mods[['element_id', 'len', 'gene_id', 'protein_id', 'start', 'end', 'feature_len', 'strand', 'target_mod','mod_rate','target_no_mod','unmodified_rate','off_target']]

#write output file
annotated_mods.to_csv(args.out+'_mod_annotations.tsv', sep='\t', index=False, header=True)

#new mod sites table with annotations
mod_sites_table = pd.DataFrame.from_dict(new_mod_sites) #*test
mod_sites_table.to_csv(args.out+'_annotated_sites.tsv', sep='\t', index=False, header=True) #*test

### summary results ###
df = annotated_mods

summary = []
for id in df.element_id.unique():

  #all mods
  filtered_df = df[df['element_id'] == id]
  total_mods = filtered_df['target_mod'].sum()
  total_len = int(filtered_df['len'].mean())
  
  #off target mods
  off_target_mods = filtered_df['off_target'].sum()
  
  #target sites not modified 
  target_no_mod = filtered_df['target_no_mod'].sum()
  
  #cds target mods
  filtered_df = df[ (df['element_id'] == id) & (df['gene_id'] != 'intergenic') ]
  cds_mods = filtered_df['target_mod'].sum()
  cds_len = filtered_df['feature_len'].sum()

  #non-coding target mods
  filtered_df = df[ (df['element_id'] == id) & (df['gene_id'] == 'intergenic') ]
  intergenic_mods = filtered_df['target_mod'].sum()
  intergenic_len = filtered_df['feature_len'].sum()

  
  #make dictionary with data 
  summary.append({'id': id,
                 'len': total_len,
                  'total_mods': total_mods,
                  'total_rate': round(total_mods/total_len*1000,3),
                 'cds_mods': cds_mods,
                 'cds_rate': round(cds_mods/cds_len*1000,3),
                 'intergenic_mods': intergenic_mods,
                  'intergenic_rate': round(intergenic_mods/intergenic_len*1000,3),
                  'off_target_mods': off_target_mods,
                  'off_target_rate': round(off_target_mods/total_len*1000,3),
                  'unmodified_targets': target_no_mod,
                  'unmodified_target_rate': round(target_no_mod/total_len*1000,3)})

#convert to df
summary_df = pd.DataFrame(data=summary)

#write output file
summary_df.to_csv(args.out+'_mod_summary.tsv', sep='\t', index=False, header=True)