#!/bin/sh

#REF=/home/firtinac/panzer/genome_remap/data/sacCer2/full/ref.fa
REF=$1
FASTQ=$2
#list of fastq suffices such as "ERR12 ERR13 ERR14" where whe have paired end reads for each
#ERRXX: ERRXX_1.fastq and ERRXX_2.fastq
BAM_SUFFIX=$3
THREADS=$4
MAXMEM=16g

for i in `echo $FASTQ`; do fname=`basename $i`; sbatch -c $THREADS --wrap="/usr/bin/time -v -p -o ${BAM_SUFFIX}_${fname}.time bwa mem -M -t $THREADS $REF ${i}\_1.fastq ${i}\_2.fastq | samtools view -h -F4 | samtools sort -m ${MAXMEM} -l0 > ${BAM_SUFFIX}_${fname}.bam"; done

