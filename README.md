This repository contains tools for analyzing 6mA modifications and building a methylome profile from Oxford Nanopore (ONT) sequencing data. Tools are built to work with the bed file output from ONT's ModKit program. 

Brief overview of workflow starting with the bed file produced by ONT's ModKit:

I. Get kmers with predicted methylated adenine in center 
  Program: get_adenine_kmers.py
  Required Inputs: 
    bed file
    ref genome
  Output: 
    windows of specified size with center adenine, fasta format
    option for control kmers to use with streme 

II. Look for enriched motifs in adenine kmers
  Program: Streme ()
  Required Inputs:
    Fasta of potentially modified adenine kmers
    Fasta of control seqs

III. Test candidate motifs for percent modified
  Programs: calculate_percent_mod.py; percent_mod_binomial_dist.py
  Required inputs:
    motif or list of motifs
    bed file
    ref genome
  Output:
    percentage of motif modified throughout genome
    binomial test statistics

IV. Map modifications to the genome
  Program: mod_mapper.py
  Inputs: 
    motif or list of motifs
    bed file
    reference genome
  Output:
    Sites table with seq id, motif, position modified, strand

V. Annotate methylated positions
  Program: annotate_methylome.py

VI. Methylation sliding window

-------------------------------------------------------

I. Synonomous codon shuffling
II. motif mapper
III. Motif enrichment avoidance analysis

Basic usage for each tool in workflow: 
