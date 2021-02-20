# AirLift

This repository contains the source code for our tool AirLift, which we describe and evaluate in our BMC Genome Biology submission, the ArXiv version (http://arxiv.org/abs/1912.08735) and the bioRxiv version (https://www.biorxiv.org/content/10.1101/2021.02.16.431517v1). 

>J.S. Kim, C. Firtina, M.B. Cavlak, D. Senol Cali, N. Hajinazar, M. Alser, C. Alkan, O. Mutlu. "AirLift: A Fast and Comprehensive Technique for Translating Alignments between Reference Genomes."

As genome sequencing tools and techniques improve, researchers are able to incrementally assemble more accurate reference genomes, which enable sensitivity in read mapping and downstream analysis such as variant calling. A more sensitive downstream analysis is critical for better understanding the health data of a genome donor. Therefore, read sets from sequenced samples should ideally be mapped to the latest available reference genome. Unfortunately, the increasingly large amount of available genomic data makes it prohibitively expensive to fully re-map each read set to its respective reference genome every time the reference is updated. There are several tools that attempt to accelerate the process of updating a read data set from one reference to another (i.e., remapping) by 1) identifying regions that appear similarly between two references and 2) updating the mapping location of reads that map to any of the identified regions in the old reference to the corresponding similar region in the new reference.  The main drawback of existing approaches is that if a read maps to a region in the old reference that does not appear similarly in the new reference, the read cannot be remapped. We find that, as a result of this drawback, a significant portion of annotations are lost when using state-of-the-art remapping tools. To address this major limitation in existing tools, we propose AirLift, a fast and comprehensive technique for moving alignments from one genome to another. AirLift reduces 1) the number of reads that need to be fully mapped from the entire read set by up to 99.99% and 2) the overall execution time to remap read sets between two reference genome versions by 19.6x, 6.6x, and 2.7x for large (human), medium (C.  elegans), and small (yeast) reference genomes, respectively. We validate our remapping results with GATK and find that AirLift provides similar rates of identifying ground truth SNPs and INDELs as fully mapping a read set.

# Running AirLift

We build AirLift on top of the following tools. Please make sure they are installed on your machine: 
* `SAMTOOLS` - A suite of programs for interacting with high-throughput sequencing data. 
* `bbmap` - A tool for finding paired-end reads within a FASTQ and splitting them into first pair, second pair, and singletons from the paired-end read data set. 
* All tools in the AirLift/dependencies directory. Install by running the install.sh script. 
* Our forked CrossMap repository at: [CrossMap](https://github.com/canfirtina/CrossMap). Be sure to include this in your PATH. 


## To run the full AirLift pipeline, please use the following command:
        
      	$ ./run_pipeline.sh [SRCFOLDER] [BINDIR] [OLDREF] [NEWREF] [SEQ_FILE_EXT] [READSIZE] [READ_BAM] [FIRST_PAIR] [SECOND_PAIR] [SAMPLE] [OUTPUT] [THREAD] [MAXMEM]

* `SRCFOLDER` - The root directory containing all the folders of different steps in AirLift. 
* `BINDIR` - The folder containing the required binaries to run AirLift. 
* `OLDREF` - The folder containing the old reference file that has already been mapped to.
* `NEWREF` - The folder containing the new reference file that AirLift will remap the read set to.
* `SEQ_FILE_EXT` - The file extension name for the reference files.
* `READSIZE` - The size of the reads that we are remapping with AirLift. 
* `READ_BAM` - The BAM files containing the reads that we are remapping with AirLift. Make sure these BAM Files are sorted and indexed (i.e., samtools index). 
* `FIRST_PAIR` - The FASTQ file containing the first pairs in the paired-end read data set. 
* `SECOND_PAIR` - The FASTQ file containing the second pairs in the paired-end read data set. 
* `SAMPLE` - The name of the sample that will be in the SAM read group header line. 
* `OUTPUT` - The folder where AirLift will place all the output files. 
* `THREAD` - The number of threads that AirLift will use to remap reads. 
* `MAXMEM` - maximum amount of memory allocated to each thread when sorting files. 



# Downstream analysis 

We call variants on AirLift output using [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) and benchmark the variant calls using [hap.py](https://github.com/Illumina/hap.py/blob/master/doc/happy.md). 

We provide an example script that runs the downstream analysis in: AirLift/src/example_airlift_script.sh

## Additional Dependencies: 

* `[MarkDuplicates-Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360036459932-MarkDuplicates-Picard-)` 


# Replicating the Results in the Paper

We provide a step-by-step guide to replicate the results reported in the paper
for Yeast. We also provide the scripts we used for C.elegans and Human genomes
in Zenodo.

Note that AirLift requires 1) an internet connection to access several
public servers (e.g., GitHub, UCSC, EBI websites), 2) python3 (and pip3), 3) 
Picard tools, 4) GATK, and 5) hap.py (https://github.com/Illumina/hap.py) 
installed on your machine.

```bash
#We will work on a directory called "test".
mkdir test; cd test;
#We clone AirLift from its GitHub page
git clone https://github.com/CMU-SAFARI/AirLift.git
#We install the required dependencies
cd AirLift/dependencies; bash install.sh; cd ../../;
#Download the Illumina paired-end read sets from the EBI website
mkdir illumina; cd illumina;
wget -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR193/003/ERR1938683/ERR1938683_1.fastq.gz | gunzip -c > ERR1938683_1.fastq
wget -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR193/003/ERR1938683/ERR1938683_2.fastq.gz | gunzip -c > ERR1938683_2.fastq
#Download the Yeast reference genome, sacCer2 and sacCer3, from the UCSC website, and create their BWA index files. We also create the Picard's dictionary files, which will be necessary for the GATK analysis.
mkdir ../sacCer2; cd ../sacCer2
wget -O - 'ftp://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/chromosomes/chr*' | gunzip -c > ref.fa
bwa index ref.fa
samtools faidx ref.fa
java -jar picard.jar CreateSequenceDictionary -R ref.fa
mkdir ../sacCer3; cd ../sacCer3;
wget -O - 'ftp://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/chromosomes/chr*' | gunzip -c > ref.fa
bwa index ref.fa
samtools faidx ref.fa
java -jar picard.jar CreateSequenceDictionary -R ref.fa
#We align the Illumina reads to the reference genomes we just downloaded
mkdir ../illumina/bwa; cd ../illumina/bwa;
/usr/bin/time -v -p -o sacCer2_ERR1938683.time bwa mem -R "@RG\tID:ERR1938683\tSM:ERR1938683\tPL:illumina\tLB:ERR1938683" -t 22 ../../sacCer2/ref.fa ../ERR1938683_1.fastq ../ERR1938683_2.fastq | samtools view -h -F4 | samtools sort -l5 -m 8G -@ 10 > sacCer2_ERR1938683.bam; samtools index sacCer2_ERR1938683.bam;
/usr/bin/time -v -p -o sacCer3_ERR1938683.time bwa mem -R "@RG\tID:ERR1938683\tSM:ERR1938683\tPL:illumina\tLB:ERR1938683" -t 22 ../../sacCer3/ref.fa ../ERR1938683_1.fastq ../ERR1938683_2.fastq | samtools view -h -F4 | samtools sort -l5 -m 8G -@ 10 > sacCer3_ERR1938683.bam; samtools index sacCer3_ERR1938683.bam;
#We prepare the sacCer3_ERR1938683.bam file for the GATK analysis that we will
#perform later (i.e., marking the PCR duplicates)
java -jar picard.jar MarkDuplicates --REMOVE_DUPLICATES -AS true -I sacCer3_ERR1938683.bam -O sacCer3_ERR1938683_rmdup.bam -M sacCer3_ERR1938683.txt
samtools index sacCer3_ERR1938683_rmdup.bam
#We create the directory to store the AirLift results and download the chain file from the UCSC website under that directory
mkdir -p ../../results/airlift; cd ../../results/airlift;
wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/liftOver/sacCer2ToSacCer3.over.chain.gz | gunzip -c > ref.chain
#We run AirLift using 32 threads, and 8GB of memory per thread used when 
#sorting (i.e., floor of 32/3 = 10 threads for sorting)
/usr/bin/time -v -p -o ./airlift.time bash ../../AirLift/src/run_pipeline.sh ../../AirLift/src/ ../../AirLift/dependencies/bin ../../sacCer2/ ../../sacCer3/ fa 150 ../../illumina/bwa/sacCer2_ERR1938683.bam ../../illumina/ERR1938683_1.fastq ../../illumina/ERR1938683_2.fastq ERR1938683 ./ 32 8G
#We prepare the AirLift-generated BAM file (i.e., ref_airlift.bam) for GATK
java -jar picard.jar CleanSam -I ref_airlift.bam -O ref_airlift_clean.bam
java -jar picard.jar MarkDuplicates --REMOVE_DUPLICATES -AS true -I ref_airlift_clean.bam -O ref_airlift_rmdup.bam -M ref_airlift.txt
rm ref_airlift_clean.bam
samtools index ref_airlift_rmdup.bam
#We run GATK HaplotypeCaller for the AirLift-generated BAM file using
#32 threads
mkdir -p ../gatk/airlift; cd ../gatk/airlift
gatk --java-options '-Xmx64G' HaplotypeCaller -R ../../../sacCer3/ref.fa -I ../../airlift/ref_airlift_rmdup.bam -O ./ref_airlift_rmdup.bam.vcf.gz --native-pair-hmm-threads 32 -RF ValidAlignmentStartReadFilter -RF ValidAlignmentEndReadFilter -OVI -ERC GVCF
bcftools view ./ref_airlift_rmdup.bam.vcf.gz | awk '$0~"^#" { print $0; next } { print $0 | "LC_ALL=C sort -k1,1 -k2,2n" }' | bcftools view -O z -o ./ref_airlift_rmdup.sorted.vcf.gz
gatk IndexFeatureFile -I ./ref_airlift_rmdup.sorted.vcf.gz
gatk GenotypeGVCFs -OVI -R ../../../sacCer3/ref.fa -V ./ref_airlift_rmdup.sorted.vcf.gz -O ./ref_airlift_rmdup.genotype.vcf.gz
#We generate the GATK results for the full mapping of reads to sacCer3
mkdir ../sacCer3; cd ../sacCer3;
gatk --java-options '-Xmx64G' HaplotypeCaller -R ../../../sacCer3/ref.fa -I ../../../illumina/bwa/sacCer3_ERR1938683_rmdup.bam -O ./sacCer3_ERR1938683_rmdup.bam.vcf.gz --native-pair-hmm-threads 32 -RF ValidAlignmentStartReadFilter -RF ValidAlignmentEndReadFilter -OVI -ERC GVCF
bcftools view ./sacCer3_ERR1938683_rmdup.bam.vcf.gz | awk '$0~"^#" { print $0; next } { print $0 | "LC_ALL=C sort -k1,1 -k2,2n" }' | bcftools view -O z -o ./sacCer3_ERR1938683_rmdup.sorted.vcf.gz
gatk IndexFeatureFile -I ./sacCer3_ERR1938683_rmdup.sorted.vcf.gz
gatk GenotypeGVCFs -OVI -R ../../../sacCer3/ref.fa -V ./sacCer3_ERR1938683_rmdup.sorted.vcf.gz -O ./sacCer3_ERR1938683_rmdup.genotype.vcf.gz
#We compare the GATK results from the AirLift-generated BAM file to the 
#baseline (i.e., full mapping) GATK results using hap.py
cd ../airlift
hap.py ../sacCer3/sacCer3_ERR1938683_rmdup.genotype.vcf.gz ./ref_airlift_rmdup.genotype.vcf.gz -r ../../../sacCer3/ref.fa -o ./ref_airlift_rmdup.genotype.vcf.gz_full_mapping --threads 32
```
