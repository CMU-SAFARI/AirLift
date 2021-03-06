#!/bin/bash

BINDIR=$1
REF=$2
FASTQ=$3
OUT_PREFIX=$4
THREADS=$5
THREAD_SORT=$6
SAMPLE=$7
MAXMEM=$8

/usr/bin/time -v -p -o ${OUT_PREFIX}.time "${BINDIR}/bwa" mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina\tLB:${SAMPLE}" -t $THREADS $REF ${FASTQ}\_1.fastq ${FASTQ}\_2.fastq | "${BINDIR}/samtools" view -h -F4 | "${BINDIR}/samtools" sort -l5 -@ ${THREAD_SORT} -m ${MAXMEM} > ${OUT_PREFIX}.bam

