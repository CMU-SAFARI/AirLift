# mrFAST

## micro-read Fast Alignment Search Tool

mrFAST is a read mapper that is designed to map short reads to reference genome with a special emphasis on the discovery of structural variation and segmental duplications. mrFAST maps short reads with respect to user defined error threshold, including indels up to 4+4 bp. This manual, describes how to choose the parameters and tune mrFAST with respect to the library settings. mrFAST is designed to find 'all'  mappings for a given set of reads, however it can return one "best" map location if the relevant parameter is invoked.

**NOTE:** mrFAST is developed for Illumina, thus requires all reads to be at the same length. For paired-end reads, lengths of mates may be different from each other, but each "side" should have a uniform length.

**Please refer to the** [wiki](https://github.com/BilkentCompGen/mrfast/wiki) **for details on how to build and use mrFAST**.

mrFAST : Micro-Read Fast Alignment Search Tool. Enhanced with FastHASH.

# Fetching and building mrFAST

Pretty simple. To fetch:

	git clone https://github.com/BilkentCompGen/mrfast.git

To build:

	cd mrfast
	make

# Usage:

	mrfast [options]

## General Options:  
	-v|--version    Current Version.  
	-h    Shows the help file.  


## Indexing Options:
	--index [file]    Generate an index from the specified fasta file.   
	--ws [int]    Set window size for indexing (default:12 max:14).  


## Searching Options:
	--search [file]    Search in the specified genome. Provide the path to the fasta file. Index file should be in the same directory.  
	--pe    Search will be done in Paired-End mode.  
	--seq [file]    Input sequences in fasta/fastq format [file]. If paired end reads are interleaved, use this option.  
	--seq1 [file]    Input sequences in fasta/fastq format [file] (First file). Use this option to indicate the first file of paired end reads.   
	--seq2 [file]    Input sequences in fasta/fastq format [file] (Second file). Use this option to indicate the second file of paired end reads.    
	-o [file]    Output of the mapped sequences. The default is "output".  
	-u [file]    Save unmapped sequences in fasta/fastq format.  
	--best    Only the best mapping from all the possible mapping is returned.  
	--seqcomp    Indicates that the input sequences are compressed (gz).  
	--outcomp    Indicates that output file should be compressed (gz).  
	-e [int]    Maximum allowed edit distance (default 4% of the read length).  
	--min [int]    Min distance allowed between a pair of end sequences.  
	--max [int]    Max distance allowed between a pair of end sequences.  
	--maxoea [int]    Max number of One End Anchored (OEA) returned for each read pair. We recommend 100 or above for NovelSeq use. Default = 100.	
	--maxdis [int]    Max number of discordant map locations returned for each read pair. We recommend 300 or above for VariationHunter use. Default = 300.  
	--crop [int]    Trim the reads to the given length.  
	--sample [string]    Sample name to be added to the SAM header (optional).  
	--rg [string]    Read group ID to be added to the SAM header (optional).  
	--lib [string]    Library name to be added to the SAM header (optional).  


## Running mrFAST via Docker

To build a Docker image:

	cd docker
	docker build . -t mrfast:latest

Your image named "mrfast" should be ready. You can run tardis using this image by

	docker run --user=$UID -v /path/to/inputs:/input -v /path/to/outputdir:/output mrfast [args]
- ```[args]``` are usual arguments you would pass to tardis executable. Be careful about mapping. You need to specify folders respective to container directory structure.
- You need to map host machine input and output directory to responding volume directories inside the container. These options are specified by '-v' argument.
- Docker works with root user by default. "--user" option saves your outputs.

Sample Docker-based command line assuming the working directory is /home/mrfast/samplerun:

	docker run --user=$UID -v /home/mrfast/samplerun/:/input -v /home/mrfast/samplerun:/output mrfast --search /input/human_g1k_v37.fasta --seq1 /input/f1.fastq --seq2 /input/f2.fastq --pe --min 0 --max 1000 -o /output/test.sam

