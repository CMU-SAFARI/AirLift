#!/bin/bash

# First clean and mark duplicates of AirLift output 
Picard-Clean_MarkDuplicates.sh <basename of AirLift output files> 

# Then run GATK Haplotyper to call variants on AirLift output 
gatk --java-options '-Xmx64G' HaplotypeCaller -OVI -R <PATH to new reference> -I <output file of Picard-Clean_MarkDuplicates.sh> -O gatk.bam.vcf.gz --native-pair-hmm-threads 16 -RF ValidAlignmentStartReadFilter -RF ValidAlignmentEndReadFilter -ERC GVCF
bcftools view gatk.bam.vcf.gz | awk '$0~"^#" { print $0; next } { print $0 | "LC_ALL=C sort -k1,1 -k2,2n" }' | bcftools view -O z -o gatk.sorted.vcf.gz 
gatk IndexFeatureFile -I gatk.sorted.vcf.gz
gatk GenotypeGVCFs -OVI -R <PATH to new reference> -V gatk.sorted.vcf.gz -O gatk.genotype.vcf.gz 


# Hap.py for analyzing GATK output of AirLift results vs. Full Mapping results against the ground truth 
hap.py <Full Mapping GATK output filenames> gatk.genotype.vcf.gz -r <PATH to new reference> -o airlift.genotype.vcf.gz_fullmapping --threads 32


