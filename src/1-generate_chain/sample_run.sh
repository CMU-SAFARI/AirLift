#!/bin/bash

#before running this command, make sure you run chain_install and install the utils folder
#somewhere. We will need this folder later to generate the chains (see below)

#this is a sample script file that assumes there are fasta files with extension .fa in both 
#folders "FROM_FA_FILES" and "TO_FA_FILES". Also it assumes that the chromosome names are 
#identical in both folders (e.g., if there is a chr1.fa in FROM_FA_FILES then there should be 
#chr.fa in TO_FA_FILES). It generates the commands to generate the chain file for each .fa file

#sample run:
#1) (optional if utils do not exist) create the utils folder by running the following command in
#folder "/path/to/folder/" bash chain_install.sh /path/to/folder/
#
#2) run this shell script assuming based on the assumptions described above
#bash sample_run.sh /path/to/script/ /path/to/binfolder /path/to/old_ref_folder/ /path/to/new_ref_folder/ "fa/fasta" /path/to/output
#The output chain files will be generated inside the folder that you run the command above (command #2)

CHAIN_SCRIPT=$1 #exact path to "chain_generate.sh"
BINDIR=$2 #the path where the bin directory folder located
FROM_FA_FILES=$3 #old reference (from)
TO_FA_FILES=$4 #new reference (to)
SEQ_FILE_EXT=$5 #what is the extension of the fasta/fastq files? possible options: fa, fasta, fastq, fq
OUTPUT=$6

for i in `echo ${FROM_FA_FILES}/*.${SEQ_FILE_EXT}`; do chr=`basename $i | sed s/.${SEQ_FILE_EXT}//`; /usr/bin/time -vpo "${OUTPUT}/$chr.time" bash ${CHAIN_SCRIPT} ${BINDIR} $i ${TO_FA_FILES}/$chr.${SEQ_FILE_EXT} "${OUTPUT}/${chr}"; mv "${OUTPUT}/${chr}/target_to_query.chain" "${OUTPUT}/${chr}.chain"; rm -rf "${OUTPUT}/${chr}"; done

