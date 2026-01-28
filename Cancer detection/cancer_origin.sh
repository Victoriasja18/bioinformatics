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

## Identify the high and moderate genes ##

# Annotate the variants calling file with the SNP sift using the reference from Gnomad
/opt/jdk-20.0.2/bin/java -jar /opt/snpEff/SnpSift.jar \
annotate ref/gnomad/af-only-gnomad.hg38.vcf.gz \
normal.freebayes.vcf > normal.fb.gnomad.vcf

# Count total number of variants 
# Grep the line that is not start with (-v) # from the annotated normal 
# Calculate each line
grep -v "^#" normal.fb.gnomad.vcf | wc -l

# Calculate including ID
# Search line by line using awk using tab ('\t') as field separator (-F) after filtering
# Select the third column for each line and if the third column is not equal to "." means have ID
grep -v "^#" normal.fb.gnomad.vcf | awk -F'\t' '$3 != "."' | wc -l

# Calculate non-ID
# Search line by line using awk using tab ('\t') as field separator (-F) after filtering
# Select the third column for each line and if the third column is equal to "." means have no ID
grep -v "^#" normal.fb.gnomad.vcf | awk -F'\t' '$3 == "."' | wc -l

# Annotate the file using snpEFF that used GRCh38 as the reference genome from the previous gnomad file 
/opt/jdk-20.0.2/bin/java -jar /opt/snpEff/snpEff.jar -v GRCh38.105 \
normal.fb.gnomad.vcf > normal.fb.snpeff.vcf

# Filter the high and moderate variants 
# Grep the line that is not start with (-v) # from the annotated normal 
grep -v "^#" normal.fb.snpeff.vcf | \
# Search line by line using awk using tab ('\t') as field separator (-F) after filtering
awk -F'\t' '{
# Look at the 8th column if they contained the annotation by ANN
  if ($8 ~ /ANN=/) {
    end=$2+length($4)-1; # Calculate the end of variants sum of the POS ($2, 2nd col) and length of ref allele - 1
    # Split the 8th cols after "ANN" into twor parts 
    # (a[1] before ANN, a[2] after ANN)
    split($8,a,"ANN=");
    # Split again into arrays anns and separates it by commas
    split(a[2],anns,",");
    # Iterate each anns 
    for(i in anns) {
    # See if the record have moderate or high
      if (anns[i] ~ /MODERATE/ || anns[i] ~ /HIGH/) {
      # Print the CHROM, POS, end, SNPid, the annotations
        print $1, $2, end, $3, anns[i];
      }
    }
  } # Print output using the field separated 
}' OFS='\t' > variant_normal.txt

less variant_normal.txt

## Check the structural variants ##

/opt/manta-1.6.0.centos6_x86_64/bin/configManta.py --bam=normal.realigned.bam \
# Using GRCh38 as reference genome 
--referenceFasta=ref/genome/Homo_sapiens_assembly38.fasta \
# Output directory for the workflow setup
# Call the regions for desired segments for call variants 
--runDir=normal_manta --callRegions=Assignment2/a2.segments.bed.gz

# Execute the workflow with 1 parallel job 
/home/s4688638/m2/a2/normal_manta/runWorkflow.py -j 1

# Look at the result 
bcftools view normal_manta/results/variants/diploidSV.vcf.gz

