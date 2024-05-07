This repository contains a workflow and tools for analyzing 6mA modifications and building a methylome profile from Oxford Nanopore (ONT) sequencing data. Tools are built to work with the bed file output from ONT's ModKit program. 

Brief overview of workflow starting with the bed file produced by ONT's ModKit:

# Motif Discovery 

I. get_adenine_kmers.py

  Retrieves the nucleotides surrounding potential modified adenines in a specified window size.
  Outputs adenine centered kmers in fasta format, along with a control set of sequences to be used to find enriched patterns 

II. Streme

  Input the positive adenine kmers and controls seqs to detect candidate methylated target motifs

III. calculate_percent_mod.py; percent_mod_binomial_dist.py

  Calculate what percentage of target sites are modified in your ONT data.
  Can perform a binomial distribuition test if necessary

-------------------------------------------------------

# Methylome Profile

I. Map modifications to the genome
  Program: mod_mapper.py
  Inputs: 
    motif or list of motifs
    bed file
    reference genome
  Output:
    Sites table with seq id, motif, position modified, strand

II. Annotate methylated positions
  Program: annotate_methylome.py

III. Methylation sliding window

-------------------------------------------------------

# Target site enrichment/avoidance analysis

I. Synonomous codon shuffling
II. motif mapper
III. Motif enrichment avoidance analysis

Basic usage for each tool in workflow: 
