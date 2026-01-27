#!/bin/bash

# De novo genome assembly

# This script performs de novo genome assembly using Shovill.
# It takes paired-end sequencing reads as input and outputs the assembled genome. 
# The SRR14457241 dataset were used for assembly from ENA database - detected as bacterial genome.
# Using illumina

# Genome assembly with shovill (using De Bruijn graph approach)
shovill --outdir ../ASSEMBLY --R1 SRR14457241.fastq.gz --R2 SRR14457241.fastq.gz

# Used this function to analysed the assembly statistics
stats.sh in=contigs.fa > contigs.fa.stats

# AB initio gene prediction with Prodigal
prodigal -a contigs.prot.faa -c -d contigs.cds.fna -f gff -g 11 -i
contigs.fa -o contigs.gff 
# -a : output protein translations
# -c : closed ends (circular)
# -d : nucleotide sequences of genes
# -f : output format (gff)
# -g : genetic code (11 for bacteria)
# -i : input contig file
# -o : output file for gene coordinates
# contigs.fa : input contig file from shovill assembly

# The assembled contigs can be further analysed for functional annotation and phylogenetic annotation.
# To do phylogenetic analysis, we required other related bacterial genomes to compare with the assembled genome
# We will used the 6 Staphylococcus family and M. caseolyticus as an outgroup for comparison.

# Using Orthofinder to indentify orthologous genes among multiple genomes
orthofinder -f ./genomes/ -t 4 -a 4