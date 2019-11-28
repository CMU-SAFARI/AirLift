#!/bin/sh

#before running this command, make sure you run chain_install and install the utils folder
#somewhere. We will need this folder later to generate the chains (see below)

#this is a sample script file that assumes there are fasta files with extension .fa in both 
#folders "FROM_FA_FILES" and "TO_FA_FILES". Also it assumes that the chromosome names are 
#identical in both folders (e.g., if there is a chr1.fa in FROM_FA_FILES then there should be 
#chr.fa in TO_FA_FILES). It generates the commands to generate the chain file for each .fa file

#sample run:
#1) create the utils folder by running the following command in folder "/path/to/folder/"
#sh chain_install.sh /path/to/folder/
#
#2) run this shell script assuming based on the assumptions described above
#sh sample_run.sh /path/to/old_ref_folder/ /path/to/new_ref_folder/ /path/to/folder/ /path/to/script/
#The output chain files will be generated inside the folder that you run the command above (command #2)

FROM_FA_FILES=$1 #old reference (from)
TO_FA_FILES=$2 #new reference (to)
UTILS=$3 #the path where the utils folder located
CHAIN_SCRIPT=$4 #exact path to "chain_generate.sh"
SEQ_FILE_EXT=$5 #what is the extension of the fasta/fastq files? possible options: fa, fasta, fastq, fq
OUTPUT=$6

for i in `echo ${FROM_FA_FILES}/*.${SEQ_FILE_EXT}`; do chr=`basename $i | sed s/.${SEQ_FILE_EXT}//`; /usr/bin/time -vpo "${OUTPUT}/$chr.time" sh ${CHAIN_SCRIPT} $i ${TO_FA_FILES}/$chr.${SEQ_FILE_EXT} ${UTILS} "${OUTPUT}/${chr}"; mv "${OUTPUT}/${chr}/target_to_query.chain" "${OUTPUT}/${chr}.chain"; rm -rf "${OUTPUT}/${chr}"; done

