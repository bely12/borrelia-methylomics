### packages needed
import pandas as pd
import operator
import collections
from collections import Counter
import Bio.Data.IUPACData as bdi
import python_codon_tables as pct
from itertools import product


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
