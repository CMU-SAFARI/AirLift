#!/bin/sh

FIRST_PAIR=$1
SECOND_PAIR=$2
BED_FILE=$3
OUTPUT=$4

/usr/bin/time -v -p -o "${OUTPUT}/extract_sequences.time" awk '{if(substr($1, length($1), 1) == 1) print substr($1, 1, length($1)-2);}' "${BED_FILE}" | seqtk subseq ${FIRST_PAIR} - > "${OUTPUT}/reads_1.fastq"
awk '{if(substr($1, length($1), 1) == 2) print substr($1, 1, length($1)-2);}' "${BED_FILE}" | seqtk subseq ${SECOND_PAIR} - > "${OUTPUT}/reads_2.fastq"
/home/firtinac/tools/release/bbmap/repair.sh in="${OUTPUT}/reads_1.fastq" in2="${OUTPUT}/reads_2.fastq" out="${OUTPUT}/fixed_reads_1.fastq" out2="${OUTPUT}/fixed_reads_2.fastq" outs="${OUTPUT}/singletons.fastq ain=t"
mv "${OUTPUT}/fixed_reads_1.fastq" "${OUTPUT}/reads_1.fastq"
mv "${OUTPUT}/fixed_reads_2.fastq" "${OUTPUT}/reads_2.fastq"

