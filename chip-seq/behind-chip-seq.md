# Behind the pipeline of ChIP-Seq

# What is ChiP-Seq?
Chromatin Immunoprecipitation Sequencing (ChIP-Seq) is a high-throughput, genome-wide method that combines antibody-based enrichment of specific DNA-protein complexes with next-generation sequencing (NGS). 

# The aims?
To identify binding sites of transcription factors, histones, and other chromatin-associated proteins

# Application
In here, we applied ChIP-Seq to identify the CTCF binding site for D. melanogaster (fruit fly). The CTCF itself is a highly conserved region with zinc-finger protein for architecture, chromatin regulation, and gene expression. 

That is why it is important to find that region to located histone modification and gene expression of the fruit fly

# What does the pipeline includes?
- Quality control using fastqc and cutadapt
- Read mapping using bowtie2 and samtools (required reference genome)
    -  samtools also to check quality mapping reads
- Peak quality check which is to prevent phantom peaks that often found in region with high sequencing depth in the ChiP-Seq analysis 
    - Using phantompeakqualtools from run_spp package
- Peak calling using macs3 
- Reproducibility analysis using idr 
- False positive filtering using bedtools 
