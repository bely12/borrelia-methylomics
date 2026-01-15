### IMPORT REQUIRED PACKAGES ###
import argparse
import random
from Bio import SeqIO
from scipy.stats import fisher_exact

### FUNCTIONS ###

# fasta reading
def load_fasta(fasta_file):
    return [str(rec.seq).upper() for rec in SeqIO.parse(fasta_file, "fasta")]

# motif matching
def motif_matches(motif, mod_pos, seq): # used in fitness scoring
    k = len(motif)
    center_index = len(seq) // 2
    start = center_index - mod_pos
    end = start + k
    if start < 0 or end > len(seq):
        return False
    window = seq[start:end]
    return all(m == b for m, b in zip(motif, window))

def motif_matches_masking(motif, mod_pos, seq): # used in motif masking 
    
    degenerate_bases = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"}}
    
    k = len(motif)
    center_index = len(seq) // 2
    start = center_index - mod_pos
    end = start + k
    if start < 0 or end > len(seq):
        return False
    window = seq[start:end]
    for m, b in zip(motif, window):
        if b not in degenerate_bases[m]:
            return False
    return True

# fitness
# def motif_fitness(motif, mod_pos, methylated, unmethylated):
#     meth_hits = sum(motif_matches(motif, mod_pos, s) for s in methylated)
#     unmeth_hits = sum(motif_matches(motif, mod_pos, s) for s in unmethylated)
#     fitness = meth_hits - unmeth_hits
#     return fitness, meth_hits, unmeth_hits

def motif_fitness(motif, mod_pos, methylated):
    k = len(motif)
    total_score = 0.0
    valid_sites = 0

    for seq in methylated:
        center_index = len(seq) // 2
        start = center_index - mod_pos
        end = start + k

        if start < 0 or end > len(seq):
            continue

        window = seq[start:end]

        matches = 0
        counted = 0

        for i, (m, b) in enumerate(zip(motif, window)):
            if i == mod_pos:
                continue  # ignore modified base
            counted += 1
            if m == b:
                matches += 1

        if counted > 0:
            total_score += matches / counted
            valid_sites += 1

    if valid_sites == 0:
        return 0.0

    return total_score / valid_sites

# mutate
def mutate_motif(motif, mod_pos, mod_base):
    bases = ["A", "C", "G", "T"]
    motif = list(motif)
    k = len(motif)
    i = random.randrange(k)
    
    if i != mod_pos:
        motif[i] = random.choice(bases)
    motif[mod_pos] = mod_base
    return "".join(motif)

# fisher exact test
def fisher_test(motif, mod_pos, methylated, unmethylated):
    pos_counts = sum(motif_matches(motif, mod_pos, seq) for seq in methylated)
    ctrl_counts = sum(motif_matches(motif, mod_pos, seq) for seq in unmethylated)
    
    table = [[pos_counts, len(methylated) - pos_counts], 
             [ctrl_counts, len(unmethylated) - ctrl_counts]]
    
    odds, pval = fisher_exact(table, alternative="greater")
    return pval, odds, pos_counts, ctrl_counts

# genetic algorithm
def evolve_motifs(methylated, unmethylated, k, mod_pos, mod_base, generations=100, pop_size=100):
    bases = ["A", "C", "G", "T"]
    population = []
    for _ in range(pop_size):
        motif = [random.choice(bases) for _ in range(k)]
        motif[mod_pos] = mod_base
        population.append("".join(motif))

    elite_size = max(1, pop_size // 10)

    for _ in range(generations):
        scored = []
        for motif in population:
            #fit, _, _ = motif_fitness(motif, mod_pos, methylated, unmethylated)
            fit = motif_fitness(motif, mod_pos, methylated)
            scored.append((motif, fit))

        scored.sort(key=lambda x: x[1], reverse=True)
        elites = [m for m, _ in scored[:elite_size]]

        newpop = elites[:]
        while len(newpop) < pop_size:
            parent = random.choice(elites)
            child = mutate_motif(parent, mod_pos, mod_base)
            newpop.append(child)

        population = newpop

    return set(population)

### ARGUMENTS ###
parser = argparse.ArgumentParser(description="Genetic algorithm for novel methylated motif discovery.")
parser.add_argument("-m", "--methylated", required=True, help='multi-fasta set of kmers centered on a known modified base, seq len needs to be odd number')
parser.add_argument("-u", "--unmethylated", required=True, help='multi-fasta set of negative control seqs, seq len needs to be odd number')
parser.add_argument("-k", "--motif_length", type=int, required=True, help='length of motif you want to search for; can search one length at a time')
parser.add_argument('-b', "--modified_base", type=str, default='A', help='the type of base modified in the motif')
parser.add_argument("-g", "--generations", type=int, default=100, help='number of generations to evolve motifs')
parser.add_argument("-n", "--population", type=int, default=100, help='number of motifs in evolving popuation')
parser.add_argument("--mask_motif", type=str, help='removes all kmers containing specified motif, useful for finding secondary, less commmon motifs')
args = parser.parse_args()

### READ POSITIVE AND NEGATIVE MULTI-FASTA SEQS ###
methylated = load_fasta(args.methylated) # pos seqs
unmethylated = load_fasta(args.unmethylated) # neg seqs

### MASK A MOTIF ###
if args.mask_motif:
    mask = args.mask_motif.upper()
    kept = []
    for s in methylated:
        found = False
        for pos in range(len(mask)):
            if mask[pos] == args.modified_base and motif_matches_masking(mask, pos, s):
                found = True
                break
        if not found:
            kept.append(s)
    methylated = kept

### RUN GENETIC ALGORITHM ###
results = []
for mod_pos in range(args.motif_length): # run alg separately for each position along length of motif 
    motifs = evolve_motifs(methylated, unmethylated, args.motif_length, mod_pos, args.modified_base, args.generations, args.population)
    
    for motif in motifs: #fisher test for likelyhood that motif is enriched in pos seqs
        p, odds, mh, uh = fisher_test(motif, mod_pos, methylated, unmethylated)
        results.append((motif, mod_pos, p, odds, mh, uh))

results.sort(key=lambda x: x[2]) # sort so smallest p val at top of table

### OUTPUT RESULTS ###
print("\nmotif\tmod_pos\tp_value\todds_ratio\tmeth_hits\tunmeth_hits") # col headers for output
for motif, pos, p, odds, mh, uh in results[:20]: # take best 20 results 
    print(f"{motif}\t{pos}\t{p:.3e}\t{odds:.2f}\t{mh}\t{uh}")

