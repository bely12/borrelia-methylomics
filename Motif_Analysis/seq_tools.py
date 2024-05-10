### packages needed
import pandas as pd
#import operator
#import heapq
#import collections
#from collections import Counter
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
    cd_table = pct.get_codons_table(139)
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
  output = []
  for seq_rec in SeqIO.parse(fasta_file, 'fasta'):
      aa_seq = seq_rec.translate()
      #if re.match(r'\*', str(aa_seq.seq)):
          #print(seq_rec.id, ": internal stop found. Skip")
          #continue
      k = 1
      for i in range(n):
          new_str = ''
          for j in range(0, len(seq_rec)-2, 3):
              codon = seq_rec[j:(j+3)]
              new_str += draw_a_syn_codon(codon)
          output.append({'id': seq_rec.id, 'seq': new_str, 'num':k})
          k+=1
  return output


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

def motif_counter(motif_list, file):
  #get all possible seqs from each motif with ambiguos nt's 
  ambig = bdi.ambiguous_dna_values
  ambig.update({'N': 'ACGT'})
  motif_dict = []
  for item in motif_list:
    all_seqs = []
    for i in product(*[ambig[j] for j in item]):
      all_seqs.append("".join(i))
    motif_dict.append({'consensus': item, 'set': all_seqs})

  seqList = list(SeqIO.parse(file, "fasta"))
  motif_counts = []
  
  for i in range(len(seqList)):
    cds = str(seqList[i].seq) #converting the seq to a string for use for regex
    rvc = str(seqList[i].seq.reverse_complement()) #making a string of the reverse complement of the seq
    
    for motif_set in motif_dict:
      for k in motif_set['set']:
        foundCDS = re.findall(pattern = r'(?=(' + k + '))', string = cds, flags = re.IGNORECASE) #searchng for a specified motif
        foundRVC = re.findall(pattern = r'(?=(' + k + '))', string = rvc, flags = re.IGNORECASE) #searchng for a specified motif

        motif_counts.append({'consensus': motif_set['consensus'],
                             'motif': k,
                             'id': seqList[i].id,
                             'seq_length': len(seqList[i].seq),
                             'fwd_match':len(foundCDS),
                             'rev_match': len(foundRVC),
                             'total_match': len(foundCDS) + len(foundRVC)})
  return motif_counts


def motif_mapper(motif_list, fasta):
  #get all possible seqs from each motif with ambiguos nt's
  ambig = bdi.ambiguous_dna_values
  ambig.update({'N': 'ACGT'})
  motif_dict = []
  for item in motif_list:
    all_seqs = []
    for i in product(*[ambig[j] for j in item]):
      all_seqs.append("".join(i))
    motif_dict.append({'consensus': item, 'set': all_seqs})

  map = []
  seqList = list(SeqIO.parse(fasta, "fasta"))
  for i in range(0, len(seqList)):
    plus_seq = str(seqList[i].seq)
    minus_seq = str(seqList[i].seq.reverse_complement())
    
    for motif_set in motif_dict:
      for k in motif_set['set']: 
        for matches in re.finditer(pattern = r'(?=(' + k + '))', string = plus_seq, flags = re.IGNORECASE): 
          map.append({'id': seqList[i].id, 'motif': motif_set['consensus'],'seq': k, 'strand': 'plus','pos': matches.start()})
        for matches in re.finditer(pattern = r'(?=(' + k + '))', string = minus_seq, flags = re.IGNORECASE): 
          map.append({'id': seqList[i].id,'motif': motif_set['consensus'], 'seq': k, 'strand': 'minus','pos': matches.start()})
  return map


def seq_replace(fasta_file_name, remove_char, replace_char, new_file_name):
  seqList = list(SeqIO.parse(fasta_file_name, "fasta"))
  new_seqs = []
  for i in range(0,len(seqList)):
    new_seqs.append(seqList[i].seq.replace(remove_char,replace_char))
    print('>',seqList[i].id,'\n',new_seqs[i],sep='', file = open(new_file_name, "a"))
