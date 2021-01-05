# AirLift

This repository contains the source code for our tool AirLift, which we describe and evaluate in our BMC Genome Biology submission and the ArXiv version (http://arxiv.org/abs/1912.08735). 

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
