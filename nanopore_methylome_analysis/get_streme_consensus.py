########## packages 
import seq_tools
import argparse
from argparse import RawTextHelpFormatter

########## ARGUMENTS
parser = argparse.ArgumentParser(
    description=''
    ' \n', formatter_class=RawTextHelpFormatter)



##### arguments for reading position weight matrix and generating consensus motif and all possible sequences
parser.add_argument('-streme_out', help='probability matrix txt file from streme output')
parser.add_argument('-motif_num', type=int, help='which  number motif from streme_out you want to use')
parser.add_argument('-get_consensus', default= False, help='takes position weight matrix and returns consensus motif wtih all possible seqs' )
parser.add_argument('-threshold',type=float, default=0.85, help='freq threshold for calling a symbol')

args = parser.parse_args()


matrix = seq_tools.get_weight_matrix(args.streme_out, args.motif_num)
matrix.index = matrix.index + 1

print(matrix)
print('\n')
consensus_motif = seq_tools.get_consensus_motif(matrix, args.threshold)
print('Consenus motif: ',consensus_motif['consensus'])

print('Possible seqs:')
for i in consensus_motif['seqs']:
  print(i)


    
