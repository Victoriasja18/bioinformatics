# Used mutect 2 to do somatic mutation caller
# Generating raw calls 
/opt/gatk-4.4.0.0/gatk Mutect2 \
-R ref/genome/Homo_sapiens_assembly38.fasta \
-I normal.realigned.bam -normal normal \
-I tumour.realigned.bam -tumor tumour \
--germline-resource ref/gnomad/af-only-gnomad.hg38.vcf.gz \
-O somatic_mut.mutect.vcf.gz \
-L Assignment2/a2.segments.bed

# Used extra germline variants from gnomad

# Generate pileups for tumour purity estimations
# Normal
/opt/gatk-4.4.0.0/gatk GetPileupSummaries \
-I normal.realigned.bam \
-V ref1/hapmap_3.3.hg38.vcf.gz \
-L ref1/hapmap_3.3.hg38.vcf.gz \
-O normal.pileups.table

# Tumour 
/opt/gatk-4.4.0.0/gatk GetPileupSummaries \
-I tumour.realigned.bam \
-V ref1/hapmap_3.3.hg38.vcf.gz \
-L ref1/hapmap_3.3.hg38.vcf.gz \
-O tumour.pileups.table

# Estimating tumour purity 
/opt/gatk-4.4.0.0/gatk CalculateContamination \
-I tumour.pileups.table \
-matched normal.pileups.table \
-O paired.contamination.table 
# NO contaminations found 

# Generate list of filtered calls 
/opt/gatk-4.4.0.0/gatk FilterMutectCalls \
-V somatic_mut.mutect.vcf.gz \
--contamination-table paired.contamination.table \
-O somatic.mutect.filtered.vcf.gz \
-R ref/Homo_sapiens_assembly38.fasta 

# Filter to the variant that only passed 
zgrep -v "^#" somatic.mutect.filtered.vcf.gz | awk -F'\t' '$7=="PASS"' >> somatic.pass.vcf

## Filtered the data usign snpEff

# Annotate using snpsift
/opt/jdk-20.0.2/bin/java -jar /opt/snpEff/SnpSift.jar \
annotate ref/gnomad/af-only-gnomad.hg38.vcf.gz \
somatic.pass.vcf > somatic.gnomad.vcf

# Annotate using snpeff
/opt/jdk-20.0.2/bin/java -jar /opt/snpEff/snpEff.jar -v GRCh38.105 \
somatic.gnomad.vcf > somatic.snpeff.vcf

# Identfy the variance that has high and moderate risk  
grep -v "^#" somatic.snpeff.vcf | \
awk -F'\t' '{
  if ($8 ~ /ANN=/) {
    end=$2+length($4)-1;
    split($8,a,"ANN=");
    split(a[2],anns,",");
    for(i in anns) {
      if (anns[i] ~ /MODERATE/ || anns[i] ~ /HIGH/) {
        print $1, $2, end, $3, anns[i];}
    }
  }
}' OFS='\t'

# Find the SV using manta
/opt/manta-1.6.0.centos6_x86_64/bin/configManta.py \
--normalBam=normal.realigned.bam \
--tumourBam=tumour.realigned.bam \
--referenceFasta=ref/genome/Homo_sapiens_assembly38.fasta \
--runDir=somatic_manta \
--callRegions=Assignment2/a2.segments.bed.gz

/home/s4688638/m2/a2/somatic_manta/runWorkflow.py -j 1

# Look at the result
bcftools view somatic_manta/results/variants/diploidSV.vcf.gz

# Extracted the chromosome, positions, reference,alternate alleles & allele frequency (AF).
bcftools query -s tumour \
  -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' \
  somatic.mutect.filtered.vcf.gz 
