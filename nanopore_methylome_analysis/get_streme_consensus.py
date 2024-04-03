########## packages 
import seq_tools
import argparse
from argparse import RawTextHelpFormatter

########## ARGUMENTS
parser = argparse.ArgumentParser(
    description='input the txt file from streme to generate consensus motifs to test for modififcaton frequency'
    ' \n', formatter_class=RawTextHelpFormatter)

parser.add_argument('-streme_out', help='probability matrix txt file from streme output')
parser.add_argument('-motif_num', type=int, help='which  number motif from streme_out you want to use')
parser.add_argument('-threshold',type=float, default=0.85, help='freq threshold for calling a symbol')
parser.add_argument('-all', const=True, nargs='?' )
args = parser.parse_args()


matrix = seq_tools.get_weight_matrix(args.streme_out, args.motif_num)
matrix.index = matrix.index + 1

print(matrix)
print('\n')
consensus_motif = seq_tools.get_consensus_motif(matrix, args.threshold)
print('Consenus motif: ',consensus_motif['consensus'])

if args.all == True:
  print('Possible seqs:')
  for i in consensus_motif['seqs']:
    print(i)


    
