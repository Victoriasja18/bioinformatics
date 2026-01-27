#!/bin/bash

# Read mapping pipeline for Lactococcus lactis Illumina data 
mkdir illum # Create new directory for illum data
# Copy the illumina data to the directory illum 
# For 0514c variant
cp Llactis_illuminaData/0514c_3_AC0617ACXX_ATCACG_L004.{1,2}.fastq 0514c.{1,2}.fastq

# For 1768int6 variant
cp illum/Llactis_illuminaData/1768int6_AC0617ACXX_TTAGGC_L004.{1,2}.fastq 1768int6.{1,2}.fastq

# Create symbolic link for the fasta file 
ln -s /opt/BINF7001/2025/Prac3_2025/Ll_MG1363.fasta Ll_MG1363.fasta

# Using Burrows-Wheeler Alignment Tool (BWA) to do alignment  
# Index the reference (reference strain MG1363)
bwa index Ll_MG1363.fasta 

# Do the alignment - using aln (generate .sai files)
mkdir aln # Create the new directory for alignment result 
# Align the 0514c read 1 with reference genome
bwa aln Ll_MG1363.fasta 0514c.1.fastq > aln/0514c.1.sai
# Align the 0514c read 2 with reference genome 
bwa aln Ll_MG1363.fasta 0514c.2.fastq > aln/0514c.2.sai 

# Align the 1768int6 read 1 with reference genome
bwa aln Ll_MG1363.fasta 1768int6.1.fastq > aln/1768int6.1.sai
# Align the 1768int6 read 2 with reference genome
bwa aln Ll_MG1363.fasta 1768int6.2.fastq > aln/1768int6.2.sai 

# Pair and map reads - using sampe (generate proper paired-end SAM files )
mkdir sam # Create new directory for mapping result 
# For each reads of 0514c based on the alignment result 
bwa sampe Ll_MG1363.fasta aln/0514c.1.sai aln/0514c.2.sai \
0514c.1.fastq \
0514c.2.fastq > sam/0514c.Ll_MG1363.sam # The output file 

# For each reads of 1768int6 based on the alignment result 
bwa sampe Ll_MG1363.fasta aln/1768int6.1.sai aln/1768int6.2.sai \
1768int6.1.fastq \
1768int6.2.fastq > sam/1768int6.Ll_MG1363.sam # The output 

# Covert and sort the sam files
mkdir bam # Create direcotry for BAM files 
# Convert the SAM to BAM for 0514c
samtools view -bS sam/0514c.Ll_MG1363.sam > bam/0514c.Ll_MG1363.bam 
# Convert the SAM to BAM for 1768int6
samtools view -bS sam/1768int6.Ll_MG1363.sam > bam/1768int6.Ll_MG1363.bam 
# Sort the BAM files For 0514c
samtools sort bam/0514c.Ll_MG1363.bam > bam/0514c.Ll_MG1363.sorted.bam 
# Sort the BAM files for 1768int6 
samtools sort bam/1768int6.Ll_MG1363.bam > bam/1768int6.Ll_MG1363.sorted.bam 

# Rename the file 
mv bam/0514c.Ll_MG1363.sorted.bam bam/0514c.Ll_MG1363.bam
mv bam/1768int6.Ll_MG1363.sorted.bam bam/1768int6.Ll_MG1363.bam
 
rm -r sam # Delete unnecessary folder

# Use read mapping 
mkdir vcf # Directory for read mapping result 
# Generate pileup of alinged reads for 0514c BAM files with reference genome
bcftools mpileup -f Ll_MG1363.fasta bam/0514c.Ll_MG1363.bam  
# Call variants from the pileup and save the output for 0154c
bcftools call --ploidy 1 -vc > vcf/0514c_vc_result.vcf 
# Generate pileup of alinged reads for 0514c BAM files againts reference genome
bcftools mpileup -f Ll_MG1363.fasta bam/1768int6.Ll_MG1363.bam
# Call variants from the pileup and save the output for 1768int6
bcftools call --ploidy 1 -vc > vcf/1768int6_vc_result.vcf
# Note: -f flag reference genome, ploidy 1 means haploid genome
# -vc use output only variant sites and create consesus variant caller 

# Compare both variant result 
# Extract specific column (2:pos, 4:ref, 5:alt, 6:qual)
cut -f2,4-6 vcf/0514c_vc_result.vcf > vcf/0514_cut.vcf # For 0514c
cut -f2,4-6 vcf/1768int6_vc_result.vcf > vcf/1768_cut.vcf # For 1768
