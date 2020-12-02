#!/bin/bash

BINDIR=$1
READ_BAM_FILE=$2
BED_FILE=$3
READSIZE=$4

awk -v bam_file="$READ_BAM_FILE" -v bin_dir="${BINDIR}/" -v read_size="${READSIZE}M" '{cmd=""bin_dir"samtools view -h "bam_file" " $1":"$2"-"$3" | "bin_dir"convert2bed --input=sam -"; while(( cmd | getline result ) > 0 ){split(result, array, "\t"); if(array[2] >= ($2-1) && array[3] <= ($3-1) && (array[5] <= 10 || array[8] != read_size) ){if(and(array[7], 64)){print array[1] "\t" array[2] "\t" array[3] "\t" array[4]".1" "\t" array[5] "\t" array[8];}else if(and(array[7], 128)){print array[1] "\t" array[2] "\t" array[3] "\t" array[4]".2" "\t" array[5] "\t" array[8];}else{print array[1] "\t" array[2] "\t" array[3] "\t" array[4] "\t" array[5] "\t" array[8];}}} close(cmd);}' ${BED_FILE} | sort -uk4,4

