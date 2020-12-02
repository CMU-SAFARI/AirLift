#!/bin/bash

BINDIR=$1
FIRST_PAIR=$2
SECOND_PAIR=$3
BED_FILE=$4
THREAD=$5
OUTPUT=$6

THREAD_SORT1=$((THREAD / 2))
THREAD_SORT2=$((THREAD - THREAD_SORT))
if [ $THREAD_SORT1 -eq 0 ]
then
	THREAD_SORT1=1
fi

awk '{if(substr($4, length($4), 1) == 1) print substr($4, 1, length($4)-2);}' "${BED_FILE}" | sort -u --parallel=${THREAD_SORT1} | "${BINDIR}/seqtk" subseq ${FIRST_PAIR} - > "${OUTPUT}/reads_1.fastq" & awk '{if(substr($4, length($4), 1) == 2) print substr($4, 1, length($4)-2);}' "${BED_FILE}" | sort -u --parallel=${THREAD_SORT2} | "${BINDIR}/seqtk" subseq ${SECOND_PAIR} - > "${OUTPUT}/reads_2.fastq"
"${BINDIR}/repair.sh" overwrite=t in="${OUTPUT}/reads_1.fastq" in2="${OUTPUT}/reads_2.fastq" out=>("${BINDIR}/rename.sh" ow=t in=stdin.fastq out="${OUTPUT}/fixed_reads_1.fastq" prefix=realigned) out2=>("${BINDIR}/rename.sh" ow=t in=stdin.fastq out="${OUTPUT}/fixed_reads_2.fastq" prefix=realigned) outs=>("${BINDIR}/rename.sh" ow=t in=stdin.fastq out="${OUTPUT}/singletons.fastq" prefix=realigned_singleton) ain=t
mv "${OUTPUT}/fixed_reads_1.fastq" "${OUTPUT}/reads_1.fastq" & mv "${OUTPUT}/fixed_reads_2.fastq" "${OUTPUT}/reads_2.fastq"

