#!/bin/bash

BINDIR=$1
READ_BAM_FILE=$2
BED_FILE=$3

awk -v bam_file="$READ_BAM_FILE" -v read_size="${READSIZE}M" '{cmd="samtools view -h "bam_file" " $1":"$2"-"$3" | convert2bed --input=sam -"; while(( cmd | getline result ) > 0 ){split(result, array, "\t"); if(array[2] >= ($2-1) && array[3] <= ($3-1) ){if(and(array[7], 64)){print array[1] "\t" array[2] "\t" array[3] "\t" array[4]".1" "\t" array[5] "\t" array[8];}else if(and(array[7], 128)){print array[1] "\t" array[2] "\t" array[3] "\t" array[4]".2" "\t" array[5] "\t" array[8];}else{print array[1] "\t" array[2] "\t" array[3] "\t" array[4] "\t" array[5] "\t" array[8];}}} close(cmd);}' ${BED_FILE} | sort -uk 4,4

