#!/bin/bash

SRCDIR=$1
BINDIR=$2
OLDREF=$3 #ref "folders" containing ref file
NEWREF=$4
SEQ_FILE_EXT=$5
READSIZE=$6
READ_BAM=$7 #make sure the bam files are sorted and indexed (i.e., samtools index)
FIRST_PAIR=$8
SECOND_PAIR=$9
SAMPLE=${10}
OUTPUT=${11}
THREAD=${12}
MAXMEM=${13} #MAXMEM per !!THREAD!! when sorting. THREAD for sorting calculated as: $THREAD/3.

#Caution: If you haven't generated the chain files yet, comment out the below line after mkdir.

mkdir -p "${OUTPUT}/"
#bash "${SRCDIR}/1-generate_chain/sample_run.sh" "${SRCDIR}/1-generate_chain/chain_generate.sh" "${BINDIR}" "${OLDREF}" "${NEWREF}" ${SEQ_FILE_EXT} "${OUTPUT}/"

for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; 
	/usr/bin/time -v -p -o "${OUTPUT}/${chr}_extracted_gaps.fa.time" python "${SRCDIR}/2-generate_gaps/extract_gaps.py" $i "${NEWREF}/${chr}.${SEQ_FILE_EXT}" ${READSIZE} "${OUTPUT}/${chr}_extracted_gaps.fa";
done

for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; 
	/usr/bin/time -v -p -o "${OUTPUT}/${chr}_kmer_gaps.fasta.time" python "${SRCDIR}/2-generate_gaps/gaps_to_fasta.py" "${OUTPUT}/${chr}_extracted_gaps.fa" ${READSIZE} "${OUTPUT}/${chr}_kmer_gaps.fasta" 10;
done

for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; 
	bash "${SRCDIR}/3-align_gaps/align_gaps.sh" "${BINDIR}" "${OLDREF}/${chr}.${SEQ_FILE_EXT}" "${OUTPUT}/${chr}_kmer_gaps.fasta" ${THREAD} "${OUTPUT}/${chr}_aligned_gaps";
done

mkdir "${OUTPUT}/bedfiles";
for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; 
	mkdir -p "${OUTPUT}/bedfiles/${chr}/";
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/extract_regions.time" bash "${SRCDIR}/4-extract_reads/extract_regions.sh" "${BINDIR}" "${OUTPUT}/${chr}_aligned_gaps.bam" "${OUTPUT}/bedfiles/${chr}/";
done

for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; 
	sed_read_scripts=$(echo ${SRCDIR}\/4-extract_reads\/extract_reads.sh | sed 's/\//\\\//g');
	sed_bindir=$(echo ${BINDIR} | sed 's/\//\\\//g');
	sed_read_bam=$(echo ${READ_BAM} | sed 's/\//\\\//g');

	ls ${OUTPUT}/bedfiles/${chr}/*.bed | sed s/$/" $sed_read_scripts $sed_bindir $sed_read_bam $READSIZE"/g | xargs -n5 sh -c 'echo "/usr/bin/time -v -p -o ${1}.reads.time bash ${2} ${3} ${4} ${1} ${5} > ${1}.reads"' sh | xargs -P ${THREAD} --replace /bin/sh -c "{}"

done

for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`;
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed.time" cat "${OUTPUT}/bedfiles/${chr}/"*.bed | "${BINDIR}/sortBed" -i - | "${BINDIR}/mergeBed" -i - > "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed";
done

for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`;
	mkdir -p "${OUTPUT}/bedfiles/${chr}/retired_bed/"; 
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed.time" python "${SRCDIR}/4-extract_reads/get_retired_regions.py" ${READSIZE} "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed" $i "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed";
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed.time" bash "${SRCDIR}/4-extract_reads/extract_reads_noprune.sh" ${BINDIR} ${READ_BAM} "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed" > "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed";
done

#constant regions:
for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`;
	mkdir -p "${OUTPUT}/bedfiles/${chr}/constant/";
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/constant/constant_regions.bed.time" python "${SRCDIR}/4-extract_reads/get_constant_regions.py" ${READSIZE} $i "${OUTPUT}/bedfiles/${chr}/constant/constant_regions.bed";
done

THREAD_SORT=$((THREAD / 3))
THREAD_MEM=$((THREAD - THREAD_SORT))
if [ $THREAD_SORT -eq 0 ]
then
	THREAD_SORT=1
fi

#realignment and merging final results
for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`;
	cat <("${BINDIR}/bwa" mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina\tLB:${SAMPLE}" "${NEWREF}/${chr}.${SEQ_FILE_EXT}" ${FIRST_PAIR} ${SECOND_PAIR} | "${BINDIR}/samtools" view -H -) <(/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/constant/constant_liftOver.time" "${BINDIR}/liftOver" -minMatch=1 <("${BINDIR}/samtools" view -h -L "${OUTPUT}/bedfiles/${chr}/constant/constant_regions.bed" ${READ_BAM} | "${BINDIR}/bamToBed" -i -) ${i} >("${SRCDIR}/5-merge/liftBedToSam" <("${BINDIR}/samtools" view -L "${OUTPUT}/bedfiles/${chr}/constant/constant_regions.bed" ${READ_BAM}) - 4 3,4 1,2) "${OUTPUT}/bedfiles/${chr}/constant/constant_unmapped.bed") | "${BINDIR}/samtools" sort -l5 -m ${MAXMEM} > "${OUTPUT}/bedfiles/${chr}/constant/constant_lifted.bam"	
	mkdir -p "${OUTPUT}/bedfiles/${chr}/realigned_reads/";
	cat <(grep ^[^#] "${OUTPUT}/bedfiles/${chr}/constant/constant_unmapped.bed") "${OUTPUT}/bedfiles/${chr}/"*.reads "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed" > "${OUTPUT}/bedfiles/${chr}/realigned_reads/realigned_reads.bed";
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/realigned_reads/extract_sequences.time" bash "${SRCDIR}/4-extract_reads/extract_sequence.sh" ${BINDIR} ${FIRST_PAIR} ${SECOND_PAIR} "${OUTPUT}/bedfiles/${chr}/realigned_reads/realigned_reads.bed" ${THREAD} "${OUTPUT}/bedfiles/${chr}/realigned_reads/";

	bash "${SRCDIR}/0-align_reads.sh" ${BINDIR} "${NEWREF}/${chr}.${SEQ_FILE_EXT}" "${OUTPUT}/bedfiles/${chr}/realigned_reads/reads" "${OUTPUT}/bedfiles/${chr}/realigned_reads/paired_reads" ${THREAD_MEM} ${THREAD_SORT} "${SAMPLE}" ${MAXMEM};
	bash "${SRCDIR}/0-align_singletons.sh" ${BINDIR} "${NEWREF}/${chr}.${SEQ_FILE_EXT}" "${OUTPUT}/bedfiles/${chr}/realigned_reads/singletons.fastq" "${OUTPUT}/bedfiles/${chr}/realigned_reads/singletons" ${THREAD_MEM} ${THREAD_SORT} "${SAMPLE}" ${MAXMEM};

	/usr/bin/time -v -p -o "${OUTPUT}/${chr}_airlift.bam.time" "${BINDIR}/samtools" merge --write-index -@ ${THREAD} --reference "${NEWREF}/${chr}.${SEQ_FILE_EXT}" -l5 -f "${OUTPUT}/${chr}_airlift.bam" "${OUTPUT}/bedfiles/${chr}/realigned_reads/paired_reads.bam" "${OUTPUT}/bedfiles/${chr}/realigned_reads/singletons.bam" "${OUTPUT}/bedfiles/${chr}/constant/constant_lifted.bam";
done

