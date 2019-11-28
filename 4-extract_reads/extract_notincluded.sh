#!/bin/bash

ORIGFSEQS=$1 #fasta or fastq file
EXISTINGREADS=$2 #assumes sorted and unique (i.e., sort | uniq)

seqtk subseq $ORIGFASTA <(comm -23 <(seqtk seq -A -C $ORIGFASTA | grep ">" | sed s/">"// | sort | uniq) $EXISTINGREADS)

