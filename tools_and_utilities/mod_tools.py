import Bio.Data.IUPACData as bdi
import python_codon_tables as pct
from itertools import product
from collections import Counter

def ambig_seq_mod(motif_dict):
  ambig = bdi.ambiguous_dna_values
  ambig.update({'N': 'ACGT'})
  all_seqs = []
  for rec in motif_dict:
    for i in product(*[ambig[j] for j in rec['motif']]):
      all_seqs.append({'seq':"".join(i), 'mod_pos': rec['mod_pos']})
  return all_seqs

def all_kmers(k):
  '''
  required for top_kmers() function
  '''
  nucleotides = ['A', 'T', 'C', 'G']
  if k < 1:
    return []
  if k == 1:
    return nucleotides
  sub_sequences = all_kmers(k - 1)
  sequences = []
  for sequence in sub_sequences:
    for nucleotide in nucleotides:
      sequences.append(nucleotide + sequence)
  return sequences

def Nmaxelements(list1, N):
    '''
    required for top_kmers() function
    '''
    final_list = []
    for i in range(0, N):
        max1 = 0
        for j in range(len(list1)):
            if list1[j] > max1:
                max1 = list1[j]
        list1.remove(max1)
        final_list.append(max1)
    return final_list

def top_kmers(seqList, mA, length_set, top_n, mod_base = 'A'):
  '''
  seqList = set of likely modified kmers with adenine in center position, multifasta file format
  mA = the position of the likely modified adenine in your input sequences ex. if they are 11mers, mA = 5
  length_set = enter as a list every kmer length you'd like to test for potential targets ex. [4,5,6] will test all kmers of length 4, 5, & 6
  top_n = top n motifs to keep/record with the greatest number of counts in your input seqs
  '''
  counts = []
  final_list = []

  #for each length target you want to search
  for k in range(length_set[0],length_set[(len(length_set)-1)]+1):

    #generate a set of all kmers of specified length
    temp = all_kmers(k)

    #keep only kmers containing at least 1 adenine
    a_kmers = []
    for kmer in temp:
      if mod_base in kmer:
        a_kmers.append(kmer)

    #one kmer at a time, iterate through candidate seqs and count the number of seqs it occurs in
    for a_kmer in a_kmers:
      count = 0
      for n in range(len(seqList)):
        seq = str(seqList[n].seq)
        #set the range within candidate seq to look in; this depends on the position of the modA and length of kmer
        if k > mA:
          sliced_seq = seq
        else:
          sliced_seq = seq[(mA+1)-k:mA+k]
        # count if match is found
        # if a_kmer in sliced_seq:
        #   count += 1
        
        new_start = mA+1 - k
        #sliced_seq = seq[new_start:mA+k]
        matched = sliced_seq.find(a_kmer) != -1
        if matched == True:
          count += 1
          og_start = new_start + sliced_seq.find(a_kmer)
          mod_pos = mA - og_start +1

      #record counts into a dictionary, then move on to the next kmer
      if count > 1:
        counts.append({'kmer_len': k, 'seq': a_kmer, 'count': count, 'mod_pos': mod_pos})

  #look at results for 1 kmer length at a time
  for length in length_set:
    filtered_counts = list(filter(lambda x: x['kmer_len'] == length, counts))

    #extract count for each kmer and put into new list
    values = []
    for i in range(len(filtered_counts)):
      values.append(filtered_counts[i]['count'])

    #identify the kmers with the n highest values and filter for kmers with those values
    top_counts = Nmaxelements(values, top_n)
    filtered_counts2 = list(filter(lambda x: x['count'] >= top_counts[len(top_counts)-1], filtered_counts))

    #print the results
    # for rec in filtered_counts2:
    #   print(rec['kmer_len'],'\t',rec['seq'],'\t',rec['count'],'\t',rec['mod_pos'])

    final_list.append(filtered_counts2)
  return final_list

def top_kmers_fast(seqList, mA, length_set, top_n, mod_base='A'):
  """
  improved version of top_kmers: only counts kmers that occur in seqList.
  """
  final_list = []

  for k in range(length_set[0], length_set[-1] + 1):
    kmer_counts = Counter()

    for record in seqList:
      seq = str(record.seq)
      if k > mA:
        sliced_seq = seq
      else:
        sliced_seq = seq[(mA+1)-k : mA+k]

      # iterate over all kmers that overlap the modified base
      for i in range(len(sliced_seq) - k + 1):
        kmer = sliced_seq[i:i+k]
        if mod_base in kmer:
          kmer_counts[kmer] += 1

    # keep only those with count > 1
    filtered_counts = [{'kmer_len': k, 'seq': seq, 'count': count} for seq, count in kmer_counts.items() if count > 1]

    # select top_n kmers
    if filtered_counts:
      top_counts = sorted(filtered_counts, key=lambda x: x['count'], reverse=True)[:top_n]
      final_list.append(top_counts)

  return final_list

def get_mod_frequency(motif_list, bed_dict, mod_pos, mod_call_thresh, min_cov = 1):
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
    #start = 6 - mod_pos # used for 11mers
    start = 16 - mod_pos # use for 31mers
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
