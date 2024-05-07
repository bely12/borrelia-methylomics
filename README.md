This repository contains a workflow and tools for analyzing 6mA modifications and building a methylome profile from Oxford Nanopore (ONT) sequencing data. Tools are built to work with the bed file output from ONT's ModKit program. 

Brief overview of workflow starting with the bed file produced by ONT's ModKit:

**Motif Discovery**

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

III. Check candidate motifs for percent modified in the genome
  Programs: calculate_percent_mod.py; percent_mod_binomial_dist.py
  Required inputs:
    motif or list of motifs
    bed file
    ref genome
  Output:
    percentage of motif modified throughout genome
    binomial test statistics

-------------------------------------------------------

Methylome Profile

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

Target site enrichment/adoidance analysis

I. Synonomous codon shuffling
II. motif mapper
III. Motif enrichment avoidance analysis

Basic usage for each tool in workflow: 
