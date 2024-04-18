import Bio.Data.IUPACData as bdi
import python_codon_tables as pct
from itertools import product

def ambig_seq_mod(motif_dict):
  ambig = bdi.ambiguous_dna_values
  ambig.update({'N': 'ACGT'})
  all_seqs = []
  for rec in motif_dict:
    for i in product(*[ambig[j] for j in rec['motif']]):
      all_seqs.append({'seq':"".join(i), 'mod_pos': rec['mod_pos']})
  return all_seqs


def get_mod_frequncy(motif_list, bed_dict, mod_pos, mod_call_thresh, min_cov = 1):
  '''
  inputs:
  motif_list is list of all possible sequences (returned in get_consesus() and dict['seqs'])
  bed_dict is the bed df converted to dictionary
  mod_pos is the predicted modified position (use its index) in the seqs
  min_cov is minimum coverage to use in search
  mod_call_threshold is the mod frequency from bed file for a position, decides if it is mod or not
  returns dictionary: {motif: n, total occurences: n, total modified: n, percent mod: n}
  '''
  results = []
  #interate through motif_list
  for motif in motif_list:
    found = 0
    modified = 0
    start = 6 - mod_pos
    end = start + len(motif)
    for i in range(len(bed_dict)):
      #if re.search(motif, bed_dict[i]['kmer'][start:end]) != None and bed_dict[i]['cov'] >= min_cov:
      if ( motif == bed_dict[i]['kmer'][start:end] ) and ( bed_dict[i]['cov'] >= min_cov ):
        found += 1
        if bed_dict[i]['mod_freq'] >= mod_call_thresh:
          modified += 1
    #calculate avg modified frequncy for motif and record in results dictionary
    if modified >= 1:
      results.append({'motif': motif, 'occurences': found, 'modified': modified, 'percent_mod': round(modified/found, 3)})
    else:
      results.append({'motif': motif, 'occurences': found, 'modified': modified, 'percent_mod': 0})
  return results


def mod_mapper(motif_dict, bed_dict, mod_call_thresh = 50.0, min_cov = 10):
  results = [] #stores all results
  used = [] #stores index from bed_dict recs that are target methylation sites 
  
  # get target methylation sites
  for rec in motif_dict:
    start = 6 - rec['mod_pos']
    end = start + len(rec['seq'])

    for i in range(len(bed_dict)):
      if ( bed_dict[i]['kmer'][start:end] == rec['seq'] ) and ( bed_dict[i]['cov'] >= min_cov ): 
        if bed_dict[i]['mod_freq'] >= mod_call_thresh:
          results.append({'motif': bed_dict[i]['kmer'][start:end], 'plasmid_id': bed_dict[i]['chromosome'], 'position': bed_dict[i]['start'], 'strand': bed_dict[i]['strand'], 'status': 'modified'})
          used.append(i)
        else:
          results.append({'motif': bed_dict[i]['kmer'][start:end], 'plasmid_id': bed_dict[i]['chromosome'], 'position': bed_dict[i]['start'], 'strand': bed_dict[i]['strand'], 'status': 'unmodified'})
          used.append(i)
  
  # get off target methylation sites
  for j in range(len(bed_dict)):
    if j in used:
      continue
    else: 
      if bed_dict[j]['mod_freq'] >= mod_call_thresh and bed_dict[j]['cov'] >= min_cov:
        results.append({'motif': 'off_target', 'plasmid_id': bed_dict[j]['chromosome'], 'position': bed_dict[j]['start'], 'strand': bed_dict[j]['strand'], 'status': 'modified'})

  return results


def mod_mapper_single(motif_list, bed_dict, mod_pos, mod_call_thresh = 50.0, min_cov = 10):
  results = []
  start = 6 - mod_pos
  end = start + len(motif_list[0])
  for i in range(len(bed_dict)):
    if ( bed_dict[i]['kmer'][start:end] in motif_list ) and ( bed_dict[i]['cov'] >= min_cov ): 
      if bed_dict[i]['mod_freq'] >= mod_call_thresh:
        results.append({'motif': bed_dict[i]['kmer'][start:end], 'plasmid_id': bed_dict[i]['chromosome'], 'position': bed_dict[i]['start'], 'strand': bed_dict[i]['strand'], 'status': 'modified'})
      else:
        results.append({'motif': bed_dict[i]['kmer'][start:end], 'plasmid_id': bed_dict[i]['chromosome'], 'position': bed_dict[i]['start'], 'strand': bed_dict[i]['strand'], 'status': 'unmodified'})
    else: 
      if ( bed_dict[i]['mod_freq'] >= mod_call_thresh ) and ( bed_dict[i]['cov'] >= min_cov ):
        if bed_dict[i]['kmer'][start:end] not in motif_list: 
          results.append({'motif': 'off_target', 'plasmid_id': bed_dict[i]['chromosome'], 'position': bed_dict[i]['start'], 'strand': bed_dict[i]['strand'], 'status': 'modified'})
  return results
