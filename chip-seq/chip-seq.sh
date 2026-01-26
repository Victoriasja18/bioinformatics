#!/bin/bash

# This pipeline aimed to find the CTCF binding site for D. melanogaster
# All of the analysis is done in the Linux command line interface inside the virtual environment created using mamba
# The data utilized in this analysis is ChIP-Seq data with 2 replicates and 1 control (input)
# The data can be downloaded from NCBI GEO database with accession number GSE47248
# Control: GSM1156583 (input)
# Treatment replicate 1: GSM1156584 (CTCF rep1)
# Treatment replicate 2: GSM1156585 (CTCF rep2) 

#=============== ENVIRONMENT SETUP===========
# We only need to do this once
# Please skip if the environment has been set up before

# Download and install Miniforge
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
# -> wget is used to download the package
bash Miniforge3-$(uname)-$(uname -m).sh
# -> bash is used to run the script, in this case, the Miniforge installer

# Create and activate a virtual environment using mamba (a more efficient version of conda)
# Replace the your_env_name with your desired name
mamba create -n [analysis] python=3.10 #this environment utilized python
# Activate to enter the ENVIRONMENT
mamba activate [analysis]

#============ Package installation =========
# We only need to do this once every creation of new ENVIRONMENT
# You can skip this step if you have installed it before
mamba install bioconda::fastqc					# To check the quality of fastq
mamba install bioconda::cutadapt				# To trim the seq; in this case, adapter
mamba install bioconda::bowtie2					# To align the sequence to the reference genome
mamba install bioconda::phantompeakqualtools	# To measure the peak
mamba install bioconda::macs3					# To differentiate the noise and the intended peak
mamba install bioconda::bedtools				# To remove false positive signals

#=========== Assess Reads Quality===========
# This step is used to check the reads quality of our sequence (FASTQ file)
# fastqc command is used to check the quality of the reads and return the html report
# Input would be fastqc to 3 FASTQ data (control, treatment_rep1, treatment_rep2)
mkdir -p fastqc_results/raw  					# Make the output folder
fastqc -t 3 -o $HOME/fastqc_results/raw $HOME/data/control.fastq.gz \
	$HOME/data/treatment_rep1.fastq.gz $HOME/data/treatment_rep2.fastq.gz

# The FastQC results shows a good data quality, except that all of the data has Overrepresented Sequence in its read. 
# In this case, this Overrepresented sequence is the adapter sequence TrueSeq, as indicated by FastQC report, with sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# To remove the adapter, we can use cutadapt or trimmomatic. In this workshop, we use cutadapt. 
# As default parameter, we set the number of thread (-j) 8, and minimum reads after trim (-m) 20
# -a parameter is used to define the adapter to be cut, -o parameter is used to define the output folder
# The input would be the sequence that we want trim with 3 fastq.gz file 

# Create directory
mkdir trimmed 	

cutadapt -j 8 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -o $HOME/trimmed/control_cutadapt.fastq.gz $HOME/data/.fastq.gz
cutadapt -j 8 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -o $HOME/trimmed/control_cutadapt.fastq.gz $HOME/data/treatment_rep1.fastq.gz
cutadapt -j 8 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -o $HOME/trimmed/control_cutadapt.fastq.gz $HOME/data/treatment_rep2.fastq.gz
# After we trimmed, we have to check the quality again. So here, we run the fastqc for all of them
mkdir -p fastqc_result/raw  					# Make the output folder
fastqc -t 3 -o $HOME/fastqc_result/raw $HOME/data/control_cutadapt.fastq.gz \
	$HOME/data/treatment_rep1_cutadapt.fastq.gz $HOME/data/treatment_rep2_cutadapt.fastq.gz
# Here, we can see that no overrepresented sequences is detected anymore
# After the data quality is good, we can continue to read mapping

#=============Read Mapping===========
# Read mapping is done to align the reads to the reference genome, so that we can understand which region these read represent.
# As such, we have to download the reference genome database first. In this workshop, the reference is already provided in $HOME/reference/ .
# Bowtie2 is the function to map the reads to genome reference. 
# In this workshop bowtie2 command, 
	# -x specifies the genome index for bowtie2,
	# -U specifies the input reads to be mapped
	# 2> redirect the log file of the mapping
# Samtools is used to sort the alignment maps and outputing the BAM format files.
# In this workshop samtools command,
	# samtools view -@ 2 -F 4 -q 1 -bS <(bowtie2 ...) is the command to convert the alignment to BAM format, so it use smaller memory
	# -F 4 is used to excludes unmapped reads (flag 4)
	# -q 1 is used to filter out non-unique mappings (MAPQ < 1)
# Bowtie and Samtools can be run using below command

# Create directory
mkdir -p mapped/bowtie2							

# Read mapping 
samtools sort -@ 2 -o $HOME/mapped/bowtie2/control_cutadapt_sorted.bam \
  <(samtools view -@ 2 -F 4 -q 1 -bS \
    <(bowtie2 -p 3 -x $HOME/reference/bowtie2_index/genome -U $HOME/trimmed/control_cutadapt.fastq.gz \
      2>$HOME/mapped/bowtie2/control.log))		# Map control reads
samtools sort -@ 2 -o $HOME/mapped/bowtie2/treatment_rep1_cutadapt_sorted.bam \
  <(samtools view -@ 2 -F 4 -q 1 -bS \
    <(bowtie2 -p 3 -x $HOME/reference/bowtie2_index/genome -U $HOME/trimmed/treatment_rep1_cutadapt.fastq.gz \
      2>$HOME/mapped/bowtie2/treatment_rep1.log)) # Map treatment_rep1 reads  
samtools sort -@ 2 -o $HOME/mapped/bowtie2/treatment_rep2_cutadapt_sorted.bam \
  <(samtools view -@ 2 -F 4 -q 1 -bS \
    <(bowtie2 -p 3 -x $HOME/reference/bowtie2_index/genome -U $HOME/trimmed/treatment_rep2_cutadapt.fastq.gz \
      2>$HOME/mapped/bowtie2/treatment_rep2.log)) 

# We can check the mapping quality using samtools flag format
samtools flagstat $HOME/mapped/bowtie2/treatment_rep1_cutadapt_sorted.bam
samtools flagstat $HOME/mapped/bowtie2/treatment_rep1_cutadapt_sorted.bam
samtools flagstat $HOME/mapped/bowtie2/treatment_rep1_cutadapt_sorted.bam

#============ Peak Quality Check =========
# To prevent phantom peaks that often found in region with high sequencing depth in the ChiP-Seq analysis
# We utilized the Phantompeakqualtools (we used package --run_spp.R) to calculates the cross-correlations scores.
# This determines if the sequences suitable for downstream analysis 
# A higher cross-correlation score indicates better ChIP enrichment and lower background noise.

# Make directory
mkdir cross_correlation

#Control 
run_spp.R -p=3 -c=$HOME/mapped/bowtie2/control_cutadapt_sorted.bam \
  -savp=$HOME/cross_correlation/control.pdf \
  -out=$HOME/cross_correlation/control.txt

#Treatment 1 
run_spp.R -p=3 -c=$HOME/mapped/bowtie2/treatment_rep1_cutadapt_sorted.bam \
  -savp=$HOME/cross_correlation/treatment_rep1.pdf \
  -out=$HOME/cross_correlation/treatment_rep1.txt

#Treatment 2
run_spp.R -p=3 -c=$HOME/mapped/bowtie2/treatment_rep2_cutadapt_sorted.bam \
  -savp=$HOME/cross_correlation/treatment_rep2.pdf \
  -out=$HOME/cross_correlation/treatment_rep2.txt

  # Input (utilized -c): the cutadapt file with .bam format 
    # We did this three times (control, treatment 1, treatment 2)
  # There are 2 possibles output that can be seen: 
  #  -savp: aimed to have the plot in .pdf format (to show the cross-correlaiton peaks)
  #  -out: output file in .txt format
      # The .txt file shows NSC (normalised strand correlations), where NSC >1.05 determines good-quality ChiP-seq data
      # It also shows RSC (relative strand correlation), where RSC > 0.8 indicated good-quality ChiP-seq 
      #   low RSC < 0.5 indicates high background noise or poor enrichment 

#============ Peak Calling ========= 
# This process is to identify the location where the protein binds to genome by analysing it's enrich region (peak) 
# In here we aimed to detect CTCF-binding sites in the D. melanogaster genome
# We utilized the package macs3 and compared the replicates with control 
  # There are 3 analysis:
    # We call both peaks from both replicates 1 and 2 against control (input both data for each replicates)
    # We call each peaks from replicates 1 and 2 against control  

# Make directory
mkdir peak_calling

# Both peaks (accounted replicate 1 and 2)
macs3 callpeak -t $HOME/mapped/bowtie2/treatment_rep1_cutadapt_sorted.bam \
               $HOME/mapped/bowtie2/treatment_rep2_cutadapt_sorted.bam \
               -c $HOME/mapped/bowtie2/control_cutadapt_sorted.bam \
               -f BAM -g dm -n both -B --outdir $HOME/peak_calling

# Peak for treatment 1
macs3 callpeak -t $HOME/mapped/bowtie2/treatment_rep1_cutadapt_sorted.bam \
               -c $HOME/mapped/bowtie2/control_cutadapt_sorted.bam \
               -f BAM -g dm -n treatment_rep1 -B --outdir $HOME/peak_calling

# Peak for treatment 2
macs3 callpeak -t $HOME/mapped/bowtie2/treatment_rep2_cutadapt_sorted.bam \
               -c $HOME/mapped/bowtie2/control_cutadapt_sorted.bam \
               -f BAM -g dm -n treatment_rep2 -B --outdir $HOME/peak_calling
  
  
  #Flags meaning:
  # -t treatment file - depends on the input
  # -c control/input file - control input
  # -f input file format (.bam)
  # -g effective genome size (dm for drosophile melanogaster) 
  # -n prefix for output name 
  # -B generate signal tracks 
  # --outdir directory for the output 

  # Format for output generated: 
    # model.r
    # peaks.narrowPeak
    # _control_lambda.bdg
    # _summits.bed
    # _pileup.bdg

#============ Reproducibility Analysis =========
# Aimed to determine if the peak are realible that detected by consistency across biological replicates
# This prevent a false positive result for the peaks
# We utilized the IDR (irreproducible discovery rate): stats method that ranks peaks based on their reproducibility
# We identify the region that have high-confidence binding sites and filter it (IDR < 0.05)

# Create new environment for reproducibility and downloading the package
mamba create -n idr_run bioconda::idr
mamba activate idr_run

# Running IDR
idr --samples $HOME/peak_calling/treatment_rep1_peaks.narrowPeak \
             $HOME/peak_calling/treatment_rep2_peaks.narrowPeak \
    --peak-list $HOME/peak_calling/both_peaks.narrowPeak \
    --input-file-type narrowPeak \
    -o $HOME/peak_calling/idrValues \
    --plot --use-best-multisummit-IDR \
    --log-output-file $HOME/peak_calling/idr.log

  #--samples the samples flag: 2 replicates peak files
  #--peak-list: pooled peak file - ranking reference
  #--input-file-type: specified the file format
  #-o output prefixes
  #--plot --use-best-multisummit-IDR : used best (lowest) IDR among multiple summits between same peak 
  #--log-output-file: path to a log file (format .log)

# Filter unrealiable IDR with threshold 0.05
awk -F'\t' '$12 > -log(0.05)/log(10)' $HOME/peak_calling/idrValues | cut -f1-10 > $HOME/peak_calling/idr.narrowPeak
  #awk -F'\t': tab as the delimiter
  #'$12 > -log(0.05)/log(10)': keep column 12 where -log10(IDR) > threshold (0.05) 
  #cut -f1-10: extract the first 10 columns that match the format 
  # > output idr.narrowpeak

#============ False positive filtering =========
# Removal of false positive signals that might be form from certain genomic regions
# We utilized bedtools 

#Download the encode blacklist (see which genomic region is giving false positive result) 
wget https://raw.githubusercontent.com/Boyle-Lab/Blacklist/61a04d2c5e49341d76735d485c61f0d1177d08a8/lists/dm6-blacklist.v2.bed.gz

#Filtering
# Keep genomic region that is not intersect with the blacklist
bedtools intersect -v -a $HOME/peak_calling/idr.narrowPeak -b dm6-blacklist.v2.bed.gz > $HOME/peak_calling/remove_blacklist
  #intersect compare 2 genomes datasets interval 
    # Input would be the the genome from our analysis (.narrowpeak) with the blacklist region (the blacklist region for D. melanogaster)
  #v inverse match returns where:
    #that returnt the entries from -a
    #that do not overlap from -b
# Our output files have reproducible peaks (IDR < 0.05) and excluding blacklist genomic region with .narrowPeak format 

# The data from the genome annotation required to be downloaded from the genome browser
# Moving data
scp ~/Downloads/dm6_tss.bed binf03@10.139.1.132:~/dm6_tss.bed

# Get CTCF binding site 
python $HOME/scripts/getclosestgene.py $HOME/peak_calling/remove_blacklist.narrowPeak $HOME/dm6_tss.bed --dist 100

#Input the narrowpeak file that already filtered
#The output would be tss_gene.bed (in .bed format)

#Cleaning the cache
mamba clean -all