#!/bin/bash
# Cancer detection and origin identification pipeline 
# All the code are run in virtual machine

# Create the symbolic limnk
ln -s /opt/BINF7001/2025/Assignment2 # For the assignment
# What is inside 
ln -s /opt/BINF7001/2025/Prac4_2025/ref . # First reference
ln -s /opt/BINF7001/2025/Prac5_2025/ref ref1 # Second reference

# We will be using the reference genome (hg38)
# Also we compared the tumor sample with normal sample

### ALIGNMENT ###
# Align the normal sample using BWA-MEM
# Add read group tags (-R) to the group
bwa mem -R "@RG\tID:BINF7001\tSM:normal\tPL:ILLUMINA" \
# Align this to the reference genome (using hg38)
ref/genome/Homo_sapiens_assembly38.fasta \
# Provide the data for the paired-end sequencing 
Assignment2/normal_1.fq.gz Assignment2/normal_2.fq.gz \
# Convert the sam to bam format
| samtools view -b - \
# Sorting the bam files and then save it into new name file
| samtools sort - > normal.sorted.bam
# View the result of alignments 
samtools flagstat normal.sorted.bam

# Align the tumour sample using BWA-MEM
# Add read group tags (-R) to the group
bwa mem -R "@RG\tID:BINF7001\tSM:patient\tPL:ILLUMINA" \
# Align this to the reference genome (using hg38)
ref/genome/Homo_sapiens_assembly38.fasta \
# Provide the data for the paired-end sequencing 
Assignment2/tumour_1.fq.gz Assignment2/tumour_2.fq.gz \
# Convert the sam to bam format 
| samtools view -b - \
# Sorting the bam files and then save it into new name file
| samtools sort - > tumour.sorted.bam
# View the result of alignments 
samtools flagstat tumours.sorted.bam

## Mark the PCR Duplicates ## 
# You can use other tools such as GATK
# For the normal sample
# Utilizing the apps piccards from java to MarkDuplicates 
/opt/java7/bin/java -jar /opt/picard/picard.jar MarkDuplicates \
# Using the interval (only focused on the specific genome)
I= normal.sorted.bam \
# Get the output mark duplicates file 
O= normal.markdups.bam \
# Get the metrics result for the mark duplicates 
M= normal.markdups.metrics
# Index the normal markduplicates
samtools index normal.markdups.bam

# For the tumour
/opt/java7/bin/java -jar /opt/picard/picard.jar MarkDuplicates \
I= tumour.sorted.bam \
O= tumour.markdups.bam \
M= tumour.markdups.metrics
samtools index tumour.markdups.bam

## Using GATK for indel realignment ##
# Normal 
# Find suspected regions around the indels 
/opt/java7/bin/java -jar /opt/GATK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-L Assignment2/a2.segments.bed \
-R ref/genome/Homo_sapiens_assembly38.fasta \
-I normal.markdups.bam \
-o normal.intervals

# Realigned indels
/opt/java7/bin/java -jar /opt/GATK/GenomeAnalysisTK.jar \
-T IndelRealigner \
-L Assignment2/a2.segments.bed \
-R ref/genome/Homo_sapiens_assembly38.fasta \
-I normal.markdups.bam \
-o normal.realigned.bam \
-targetIntervals normal.intervals \
--disable_bam_indexing 
samtools index normal.realigned.bam

# Tumour 
# Find suspected regions around the indels 
/opt/java7/bin/java -jar /opt/GATK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-L Assignment2/a2.segments.bed \
-R ref/genome/Homo_sapiens_assembly38.fasta \
-I tumour.markdups.bam \
-o tumour.intervals

# Realigned indels
/opt/java7/bin/java -jar /opt/GATK/GenomeAnalysisTK.jar \
-T IndelRealigner \
-L Assignment2/a2.segments.bed \
-R ref/genome/Homo_sapiens_assembly38.fasta \
-I tumour.markdups.bam \
-o tumour.realigned.bam \
-targetIntervals tumour.intervals \
--disable_bam_indexing 
samtools index tumour.realigned.bam

## Variant Calling ##
# Using freebayes
# Normal
freebayes -F 0.2 --min-repeat-entropy 0 -f ref/genome/Homo_sapiens_assembly38.fasta \
normal.realigned.bam > normal.freebayes.vcf

#Tumour
freebayes -F 0.2 --min-repeat-entropy 0 -f ref/genome/Homo_sapiens_assembly38.fasta \
tumour.realigned.bam > tumour.freebayes.vcf