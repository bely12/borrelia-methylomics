This repository contains a workflow and tools for analyzing 6mA modifications and building a methylome profile from Oxford Nanopore (ONT) sequencing data. Tools are built to work with the bed file output from ONT's ModKit program. 

For a detailed tutorial for using these tools (with code), visit https://bely12.github.io/

Brief overview of workflow starting with the bed file produced by ONT's ModKit:

# Motif Discovery 

This module scans methylation calls from ONT's ModKit and puts together DNA consensus motifs targeted for methylation in your genome, and calculates genome wide methylation frequencies for the motifs. Analysis can be expanded with statistical tests on enrichment and avoidance of the motifs in your reference genome. 

-------------------------------------------------------

# Methylome Profiler

Maps methylated motifs to a reference sequence, annotates methylated positions with overlapping genomic features from a gbff file, and performs a sliding window analysis to calculate methylation frequency across intervals within a genome.
