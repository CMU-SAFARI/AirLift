#!/bin/bash

#takes the gaps and old reference file and aligns the gaps to the reference with 2x edit distance of
#the default edit distance value (i.e., 8%)
#sample run:
#sh align_gaps.sh /path/to/index/oldref.fa /path/to/gaps/gaps.fa 48 /path/to/output/prefix

BINDIR=$1
INDEX=$2 #prefix of the index files
GAPS=$3
THREAD=$4 #number of threads allocated for this job
OUTPUT_PREFIX=$5

/usr/bin/time -v -p -o "${OUTPUT_PREFIX}_aln.time" "bwa" aln -n 0.08 -t "${THREAD}" "${INDEX}" "${GAPS}" > "${OUTPUT_PREFIX}.sai"
/usr/bin/time -v -p -o "${OUTPUT_PREFIX}_samse.time" "bwa" samse "${INDEX}" "${OUTPUT_PREFIX}.sai" "${GAPS}" | "samtools" view -h -F4 | "samtools" sort -l5 -m32G > "${OUTPUT_PREFIX}.bam"

