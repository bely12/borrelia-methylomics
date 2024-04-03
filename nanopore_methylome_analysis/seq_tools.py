#TO Do 
#for get_consensus function, make row selection for fetching matrix better 
#can probably make random_seq and custom_random_seq into one function

import pandas as pd
import operator
import heapq
import collections
from collections import Counter
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.Data.IUPACData as bdi
import python_codon_tables as pct
from itertools import product
import re
import random


def ambig_seq(seq):
  '''
  input a sequence with ambiguous nt's
  returns a list of all possible sequences
  '''
  ambig = bdi.ambiguous_dna_values
  ambig.update({'N': 'ACGT'})
  all_seqs = []
  for i in product(*[ambig[j] for j in seq]):
    all_seqs.append("".join(i))
  return all_seqs

def get_weight_matrix(streme_txt, motif_num): #for use with streme .txt output
  '''
  input the streme output txt file with the position weight matrices; motif_num is the number of the motif you want
  returns a dataframe of the weight matrix of the specified motif
  '''
  x = pd.read_fwf(streme_txt, header=None) #upload streme prob matrix txt file
  index_dict = []
  m = 1 #for labeling motifs in the index_dict (first motif = 1 etc.)
  cycle = 0 #used to skip to next matrix in file after the first one is mapped
  for i in range(len(x[0])-1): #has to be -1 to avoid being out of range since for j loop looks at next line from i
    i = cycle # used to avoid reading every line, but rather find the first matrix, log it, then move to the next one
    if x[0][i].startswith('0') or x[0][i].startswith('1') == True: #conditional to find boundaries of matrix in file
      start = i #represents the first line with the matrix table
      k = 1 #used to look at the next line after i
      for j in range(len(x[0])): #only needs to iterate the length of the motif, but it varies so just use full len of x[0] and break after matrix
        if x[0][i+k].startswith('0') == False and x[0][i+k].startswith('1') == False: #conditional as above but inverse to find end
          end = i + k #represents last line in txt file for weight matrix
          index_dict.append({'motif':m, 'start':start, 'end':end}) #store coordinates of weight matrix in txt file
          cycle = end + 1 #this ensures that the "for i" loop jumps to the position we left off at
          m +=1 # moving on to the next motif in the file for the next loop
          break #break out of "for j..." loop now that weight matrix coordinates have been logged - goes back to "for i..." loop
        else:
          k +=1 #keep the "for j..." loop going but increases k to go to next line
    else:
      cycle = i + 1 #not sure if necessary here but it works
    if cycle >= len(x[0]): #stop when end of file is reached
      break
  #build dataframe
  pos1 = index_dict[motif_num -1]['start'] # use -1 to make more user friendly arg in function (1st motif is 1 instread of 0)
  pos2 = index_dict[motif_num -1]['end']
  weight_matrix = x.iloc[pos1:pos2, 0] #extract frequencies for first motif
  weight_matrix = pd.DataFrame(weight_matrix) #convert to a datafram
  weight_matrix[['A','C','G','T']] = weight_matrix[0].str.split(' ', expand=True) #add col names for nt's, split lines into col by space
  weight_matrix = weight_matrix.drop(columns=[0]).reset_index(drop=True) #drop first col not needed anymore
  return weight_matrix

def get_consensus_motif(weight_matrix, threshold = 0.85):
  '''
  input position weight matrix as a pandas df, specify threshold for calling a base at a position
  returns a dictionary: {consesus: consensus motif, seqs: all possible sequences for consensus motif}
  * requires ambig_seq() function *
  '''
  weight_matrix = weight_matrix.to_dict('records') #convert weight matrix to dicitonary
  ambig_nt = bdi.ambiguous_dna_values #get dictionary of ambiguous bases
  ambig_nt.update({'N': 'ACGT'}) #change order of bases for N to make easier to match later (using alphabetical order)

  #algorithm
  seqs = [] #will become a list, each item is string of all possible nt's at each position
  for position in range(len(weight_matrix)): #position is index for list of dictionaries
    if float(weight_matrix[position].get(max(weight_matrix[position].items(), key=operator.itemgetter(1))[0])) >= threshold: #does nt with highest freq exceed threshold?
      seqs.append(max(weight_matrix[position].items(), key=operator.itemgetter(1))[0]) #add NT to list
    else:
      nt2 = Counter(weight_matrix[position]).most_common(2) #get 2 most freq bases
      if float(nt2[0][1]) + float(nt2[1][1]) >= threshold: #check if 2 highest freqs exceed threshold
        temp = nt2[0][0] + nt2[1][0] #combine both as a double nt
        temp2 = ''.join(sorted(temp)) #double nt to alphabetical order
        seqs.append(temp2) #add string to list
      else:
        nt3 = Counter(weight_matrix[position]).most_common(3) #same as above but for triple nt
        if float(nt3[0][1]) + float(nt3[1][1]) + float(nt3[2][1]) >= threshold:
          temp = nt3[0][0] + nt3[1][0] + nt3[2][0]
          temp2 = ''.join(sorted(temp))
          seqs.append(temp2)
        else:
          seqs.append('ACGT') #N if not single, double, or triple
  #seqs
  consensus_motif = '' #start empty string for consensus seq with ambiguous symbols
  for i in range(len(seqs)):
    for key in ambig_nt.keys():
      if ambig_nt[key] == seqs[i]:
        consensus_motif = consensus_motif + key #conver any double/triple/quadruple nt to ambiguous symbols
        break
      else:
        continue
  #print(consensus_motif)
  motif_dict = {'consensus': consensus_motif, 'seqs': ambig_seq(consensus_motif)}
  return motif_dict

def get_mod_frequncy(motif_list, bed_dict, mod_pos, mod_call_thresh, min_cov = 1): #for use with bed file from mCaller
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
  #interate through motif_list, clear mod_freq list for each motif
  for motif in motif_list:
    mod_freq = []
    found = 0
    modified = 0
    start = 6 - mod_pos
    end = start + len(motif)
    for i in range(len(bed_dict)):
      if re.search(motif, bed_dict[i]['kmer'][start:end]) != None and bed_dict[i]['cov'] >= min_cov:
        found += 1
        if bed_dict[i]['mod_freq'] >= mod_call_thresh:
          modified += 1
    #calculate avg modified frequncy for motif and record in results dictionary
    if modified >= 1:
      results.append({'motif': motif, 'occurences': found, 'modified': modified, 'percent_mod': round(modified/found, 3)})
    else:
      results.append({'motif': motif, 'occurences': found, 'modified': modified, 'percent_mod': 0})
  return results

def base_composition(seq):
  A = seq.count('a')
  T = seq.count('t')
  C = seq.count('c')
  G = seq.count('g')
  composition = A/len(seq)*100, T/len(seq)*100, C/len(seq)*100, G/len(seq)*100
  return composition

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'M': 'K', 'K': 'M'}
    #complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    #bases = list(seq.lower())
    bases = list(seq.upper())
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def reverse_complement(s):
    return complement(s[::-1])

def draw_a_syn_codon(codon):
    aa = codon.translate() # aa residue as a seqRecord
    syn_codons = cd_table[str(aa.seq)]
    probs = list(syn_codons.values())
    #print(probs, end = ";")
    probs[-1] = 1 - sum(probs[0:-1]) # force sum to 1
    #print(probs)
    new_codon = np.random.choice(list(syn_codons.keys()), size = 1, replace = False, p = probs)
    return new_codon[0]
    #print(syn_codons)

def syn_codon_shuffle(fasta_file, n):
    output = [] #added this much later, doesn't make sense that it would work without it? 
    for seq_rec in SeqIO.parse(fasta_file, 'fasta'):
        aa_seq = seq_rec.translate()
        #if re.match(r'\*', str(aa_seq.seq)):
            #print(seq_rec.id, ": internal stop found. Skip")
            #continue
        for i in range(n):
            new_str = ''
            for j in range(0, len(seq_rec)-2, 3):
                codon = seq_rec[j:(j+3)]
                new_str += draw_a_syn_codon(codon)
            output.append({'id': seq_rec.id, 'seq': new_str})
        #for item in output:
            #print('>',item['id'],'\n',item['seq'],sep='')

def random_seq(seq_len):
  bases = ['A', 'T', 'C', 'G']
  output = []
  for i in range(0,seq_len):
    x = (random.choice(bases))
    output.append(x)
  random_seq = ''.join(output)
  return random_seq 

def custom_random_seq(item_bank, weights, seq_len):
  random_seq = (random.choices(item_bank, weights = weights, k = seq_len))
  return random_seq

def weighted_seq_generator(wt_seq): #input list of nt's and wt_seq as a string
  comp = base_composition(wt_seq)
  bases = ['a', 't', 'c', 'g']
  random_seq = (random.choices(bases, weights = comp, k = len(wt_seq)))
  return random_seq

def seq_shuffler(file, n):
  seqList = list(SeqIO.parse(file, "fasta"))
  for i in range(len(seqList)):
    #k = 1
    for z in range(n):
      seq=list(seqList[i].seq)
      random.shuffle(seq)
      seq = ''.join(seq)
      #output.append({'id': seqList[i].id+'_shuffled_'+str(k), 'seq': seq}) #use for unique sequence names
      output.append({'id': seqList[i].id, 'seq': seq})
      #k = k+1

def motif_counter(motif, file):
    seqList = list(SeqIO.parse(file, "fasta"))
    for i in range(0, len(seqList)):
        cds = str(seqList[i].seq) #converting the seq to a string for use for regex
        rvc = str(seqList[i].seq.reverse_complement()) #making a string of the reverse complement of the seq
        foundCDS = re.findall(pattern = r'(?=(' + motif + '))', string = cds, flags = re.IGNORECASE) #searchng for a specified motif
        foundRVC = re.findall(pattern = r'(?=(' + motif + '))', string = rvc, flags = re.IGNORECASE) #searchng for a specified motif
        output.append({'id': seqList[i].id, 'seq_length': len(seqList[i].seq), 'motif': motif, 
                   'fwd_match':len(foundCDS), 
                   'rev_match': len(foundRVC),
                  'total_match': len(foundCDS) + len(foundRVC)})

def motif_mapper(motif, file):
    seqList = list(SeqIO.parse(file, "fasta"))
    for i in range(0, len(seqList)):
      plus_seq = str(seqList[i].seq)
      minus_seq = str(seqList[i].seq.reverse_complement()) 
      for matches in re.finditer(pattern = r'(?=(' + motif + '))', string = plus_seq, flags = re.IGNORECASE): 
        output.append({'id': seqList[i].id, 'seq_length': len(seqList[i].seq), 'motif': motif, 'strand': 'plus','pos_start': matches.start()})
      for matches in re.finditer(pattern = r'(?=(' + motif + '))', string = minus_seq, flags = re.IGNORECASE): 
        output.append({'id': seqList[i].id, 'seq_length': len(seqList[i].seq), 'motif': motif, 'strand': 'minus','pos_start': matches.start()})

def seq_replace(fasta_file_name, remove_char, replace_char, new_file_name):
  seqList = list(SeqIO.parse(fasta_file_name, "fasta"))
  new_seqs = []
  for i in range(0,len(seqList)):
    new_seqs.append(seqList[i].seq.replace(remove_char,replace_char))
    print('>',seqList[i].id,'\n',new_seqs[i],sep='', file = open(new_file_name, "a"))
