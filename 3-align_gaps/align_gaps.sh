#!/bin/sh

#takes the gaps and old reference file and aligns the gaps to the reference with 2x edit distance of
#the default edit distance value (i.e., 8%)
#sample run:
#sh align_gaps.sh /path/to/index/oldref.fa /path/to/gaps/gaps.fa 48 /path/to/output/prefix

INDEX=$1 #prefix of the index files
GAPS=$2
THREAD=$3 #number of threads allocated for this job
OUTPUT_PREFIX=$4

/usr/bin/time -v -p -o "${OUTPUT_PREFIX}_aln.time" bwa aln -n 0.08 -t "${THREAD}" "${INDEX}" "${GAPS}" > "${OUTPUT_PREFIX}.sai"
/usr/bin/time -v -p -o "${OUTPUT_PREFIX}_samse.time" bwa samse "${INDEX}" "${OUTPUT_PREFIX}.sai" "${GAPS}" | samtools view -h -F4 | samtools sort -m 16G -l0 > "${OUTPUT_PREFIX}.bam"

