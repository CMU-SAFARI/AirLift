# AirLift

This repository contains the source code for our tool AirLift, which we describe and evaluate in our RECOMB'20 submission (under review) and the ArXiv version (http://arxiv.org/abs/1912.08735). 

>J.S. Kim, C. Firtina, D. Senol Cali, M. Alser, N. Hajinazar, C. Alkan, O. Mutlu. "AirLift: A Fast and Comprehensive Technique for Translating Alignments between Reference Genomes."

As genome sequencing tools and techniques improve, researchers are able to incrementally assemble more accurate reference genomes. A more accurate reference genome enables increased accuracy in read mappings, which provides more accurate  variant information and thus health data on the donor. Therefore, read data sets from sequenced samples should ideally be mapped to the latest available reference genome. Unfortunately, the increasingly large amounts of available genomic data makes it prohibitively expensive to fully map each read data set to its respective reference genome every time the reference is updated. Several tools that attempt to reduce the procedure of updating a read data set from one reference to another (i.e., remapping) have been published. These tools identify regions of similarity across the two references and update the mapping locations of a read based on the locations of similar regions in the new reference genome. The main drawback of existing approaches is that if a read maps to a region in the old reference without similar regions in the new reference, it cannot be remapped. We find that, as a result of this drawback, a significant portion of annotations are lost when using state-of-the-art remapping tools. To address this major limitation in existing tools, we propose AirLift, a fast and comprehensive technique for moving alignments from one genome to another. AirLift can reduce 1) the number of reads that need to be mapped from the entire read set by up to 99.9\% and 2) the overall execution time to remap the reads between the two most recent reference versions by 6.94$\times$, 44.0$\times$, and 16.4$\times$ for large (human), medium (C. elegans), and small (yeast) references, respectively.

# Running AirLift

We build AirLift on top of the following tools. Please make sure they are installed on your machine: 
* `SAMTOOLS` - A suite of programs for interacting with high-throughput sequencing data. 
* `slurm` - A workload manager for job scheduling on a compute cluster. 
* `bbmap` - A tool for finding paired-end reads within a FASTQ and splitting them into first pair, second pair, and singletons from the paired-end read data set. 


## To run AirLift, please issue the following command:
        
       
      	$ ./run_pipeline.sh [SRCFOLDER] [OLDREF] [NEWREF] [SEQ_FILE_EXT] [READSIZE] [READ_BAM] [FIRST_PAIR] [SECOND_PAIR] [OUTPUT] [THREAD] 
        

* `SRCFOLDER` - The root directory containing all the folders of different steps in AirLift. 
* `OLDREF` - The folder containing the old reference file.
* `NEWREF` - The folder containing the new reference file.
* `SEQ_FILE_EXT` - The file extension name for the reference files.
* `READSIZE` - The size of the reads that we are remapping with AirLift. 
* `READ_BAM` - The BAM files containing the reads that we are remapping with AirLift. Make sure these BAM Files are indexed (i.e., samtools index). 
* `FIRST_PAIR` - The FASTQ file containing the first pairs in the paired-end read data set. 
* `SECOND_PAIR` - The FASTQ file containing the second pairs in the paired-end read data set. 
* `OUTPUT` - The folder where AirLift will place all the output files. 
* `THREAD` - The number of threads that AirLift will use to remap reads. 

## This will perform all the steps of AirLift from end to end: 

1. `generate_chain` - generate a chainfile, or a file containing the similar regions across two reference genomes. This step does not need to be taken if a chain file is provided for the pair of references. 
2. `generate_gaps` - create a FASTA file containing the gaps, or regions of the old reference genome that does not exist in the new reference genome (from the information encoded in the chain file). 
3. `align_gaps` - align k-mers from the gaps to the old reference genome. 
4. `extract_reads` - extracts reads that mapped to different types of regions in the old reference genome that we need to map to the new reference genome using hints provided by our maps. 

Note that steps 1-3 only need to be run once for a pair of reference genomes to remap any number of reads between the two reference genomes.  

