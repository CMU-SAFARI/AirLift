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

#Caution: If you haven't generated the chain files yet, uncomment the command immediately following the mkdir command (i.e., chain_generate.sh).

# Create folder where all output files will be placed 
mkdir -p "${OUTPUT}/"

# Build the chain file. Uncomment this line if you do not already have a chain file, for the pair of references you want to remap reads between. 
#bash "${SRCDIR}/1-generate_chain/sample_run.sh" "${SRCDIR}/1-generate_chain/chain_generate.sh" "${BINDIR}" "${OLDREF}" "${NEWREF}" ${SEQ_FILE_EXT} "${OUTPUT}/"

#######################################################################################
# The next set of commands use the chain file to categorize each region in the references as constant, updated, retired, or new. The reads from these regions are then separated in files according to their corresponding region's categorization. 
#######################################################################################

# Create a fasta file only containing the non-constant regions in the new reference as indicated by the chain file. 
for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; 
	/usr/bin/time -v -p -o "${OUTPUT}/${chr}_extracted_gaps.fa.time" python "${SRCDIR}/2-generate_gaps/extract_gaps.py" $i "${NEWREF}/${chr}.${SEQ_FILE_EXT}" ${READSIZE} "${OUTPUT}/${chr}_extracted_gaps.fa";
done

# Extract all overlapping read-sized tokens from the fasta file containing non-constant regions
for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; 
	/usr/bin/time -v -p -o "${OUTPUT}/${chr}_kmer_gaps.fasta.time" python "${SRCDIR}/2-generate_gaps/gaps_to_fasta.py" "${OUTPUT}/${chr}_extracted_gaps.fa" ${READSIZE} "${OUTPUT}/${chr}_kmer_gaps.fasta" 10;
done

# align read-sized tokens to the old reference using BWA 
for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; 
	bash "${SRCDIR}/3-align_gaps/align_gaps.sh" "${BINDIR}" "${OLDREF}/${chr}.${SEQ_FILE_EXT}" "${OUTPUT}/${chr}_kmer_gaps.fasta" ${THREAD} "${OUTPUT}/${chr}_aligned_gaps";
done

# Get all regions of the old reference that the read-sized tokens mapped to
mkdir "${OUTPUT}/bedfiles";
for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; 
	mkdir -p "${OUTPUT}/bedfiles/${chr}/";
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/extract_regions.time" bash "${SRCDIR}/4-extract_reads/extract_regions.sh" "${BINDIR}" "${OUTPUT}/${chr}_aligned_gaps.bam" "${OUTPUT}/bedfiles/${chr}/";
done

# Extract all the reads from the regions in the old reference that the read-sized tokens mapped to (identified in the prior step) 
for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; 
	sed_read_scripts=$(echo ${SRCDIR}\/4-extract_reads\/extract_reads.sh | sed 's/\//\\\//g');
	sed_bindir=$(echo ${BINDIR} | sed 's/\//\\\//g');
	sed_read_bam=$(echo ${READ_BAM} | sed 's/\//\\\//g');

	ls ${OUTPUT}/bedfiles/${chr}/*.bed | sed s/$/" $sed_read_scripts $sed_bindir $sed_read_bam $READSIZE"/g | xargs -n5 sh -c 'echo "/usr/bin/time -v -p -o ${1}.reads.time bash ${2} ${3} ${4} ${1} ${5} > ${1}.reads"' sh | xargs -P ${THREAD} --replace /bin/sh -c "{}"
done

for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`;
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed.time" cat "${OUTPUT}/bedfiles/${chr}/"*.bed | "${BINDIR}/sortBed" -i - | "${BINDIR}/mergeBed" -i - > "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed";
done

# identify the retired regions (i.e., regions that are not similar to any region in the new reference) in the old reference and extract the reads from those regions. 
for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`;
	mkdir -p "${OUTPUT}/bedfiles/${chr}/retired_bed/";
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed.time" python "${SRCDIR}/4-extract_reads/get_retired_regions.py" ${READSIZE} "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed" $i "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed";
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed.time" bash "${SRCDIR}/4-extract_reads/extract_reads_noprune.sh" ${BINDIR} ${READ_BAM} "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed" > "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed";
done

#######################################################################################
# remap reads from the constant regions using CrossMap. 
#######################################################################################

# Identify the constant regions based on the chain file. 
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

# Remap constant regions from old reference to new reference using CrossMap
for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`;
	"${BINDIR}/samtools" view -hb -@ ${THREAD} -L "${OUTPUT}/bedfiles/${chr}/constant/constant_regions.bed" ${READ_BAM} > "${OUTPUT}/bedfiles/${chr}/constant/crossmap_before.bam"
	"${BINDIR}/samtools" index -@ ${THREAD} "${OUTPUT}/bedfiles/${chr}/constant/crossmap_before.bam"
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/constant/crossmap_final.time" CrossMap.py bam $i "${OUTPUT}/bedfiles/${chr}/constant/crossmap_before.bam" "${OUTPUT}/bedfiles/${chr}/constant/crossmap_final.bed" - -a | "${BINDIR}/samtools" sort -@ ${THREAD_SORT} -l5 -m ${MAXMEM} > "${OUTPUT}/bedfiles/${chr}/constant/crossmap_final.bam"
	"${BINDIR}/samtools" index -@ ${THREAD} "${OUTPUT}/bedfiles/${chr}/constant/crossmap_final.bam"
	rm "${OUTPUT}/bedfiles/${chr}/constant/crossmap_before.bam" "${OUTPUT}/bedfiles/${chr}/constant/crossmap_before.bam.bai"
	#cat <("${BINDIR}/bwa" mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina\tLB:${SAMPLE}" "${NEWREF}/${chr}.${SEQ_FILE_EXT}" ${FIRST_PAIR} ${SECOND_PAIR} | "${BINDIR}/samtools" view -H -) <(/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/constant/constant_liftOver.time" "${BINDIR}/liftOver" -minMatch=1 <("${BINDIR}/samtools" view -h -L "${OUTPUT}/bedfiles/${chr}/constant/constant_regions.bed" ${READ_BAM} | "${BINDIR}/bamToBed" -i -) ${i} >("${SRCDIR}/5-merge/liftBedToSam" <("${BINDIR}/samtools" view -L "${OUTPUT}/bedfiles/${chr}/constant/constant_regions.bed" ${READ_BAM}) - 4 3,4 1,2) "${OUTPUT}/bedfiles/${chr}/constant/constant_unmapped.bed") | "${BINDIR}/samtools" sort -l5 -m ${MAXMEM} > "${OUTPUT}/bedfiles/${chr}/constant/constant_lifted.bam"	
	mkdir -p "${OUTPUT}/bedfiles/${chr}/realigned_reads/";
	#cat <(grep ^[^#] "${OUTPUT}/bedfiles/${chr}/constant/constant_unmapped.bed")
	cat "${OUTPUT}/bedfiles/${chr}/constant/crossmap_final.bed" "${OUTPUT}/bedfiles/${chr}/"*.reads "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed" > "${OUTPUT}/bedfiles/${chr}/realigned_reads/realigned_reads.bed";
	/usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/realigned_reads/extract_sequences.time" bash "${SRCDIR}/4-extract_reads/extract_sequence.sh" ${BINDIR} ${FIRST_PAIR} ${SECOND_PAIR} "${OUTPUT}/bedfiles/${chr}/realigned_reads/realigned_reads.bed" ${THREAD} "${OUTPUT}/bedfiles/${chr}/realigned_reads/";

#######################################################################################
# Align the updated and retired reads to the new reference 
#######################################################################################
	bash "${SRCDIR}/0-align_reads.sh" ${BINDIR} "${NEWREF}/${chr}.${SEQ_FILE_EXT}" "${OUTPUT}/bedfiles/${chr}/realigned_reads/reads" "${OUTPUT}/bedfiles/${chr}/realigned_reads/paired_reads" ${THREAD_MEM} ${THREAD_SORT} "${SAMPLE}" ${MAXMEM};
	bash "${SRCDIR}/0-align_singletons.sh" ${BINDIR} "${NEWREF}/${chr}.${SEQ_FILE_EXT}" "${OUTPUT}/bedfiles/${chr}/realigned_reads/singletons.fastq" "${OUTPUT}/bedfiles/${chr}/realigned_reads/singletons" ${THREAD_MEM} ${THREAD_SORT} "${SAMPLE}" ${MAXMEM};

#######################################################################################
# Merge all the files containing remapped reads for a final AirLift output file 
#######################################################################################
	cat <("${BINDIR}/samtools" view -H "${OUTPUT}/bedfiles/${chr}/constant/crossmap_final.bam") <("${BINDIR}/samtools" merge -@ ${THREAD_SORT} - "${OUTPUT}/bedfiles/${chr}/realigned_reads/paired_reads.bam" "${OUTPUT}/bedfiles/${chr}/realigned_reads/singletons.bam" | "${BINDIR}/samtools" view) | "${BINDIR}/samtools" sort -l5 -m ${MAXMEM} -@ ${THREAD_SORT} > "${OUTPUT}/bedfiles/${chr}/realigned_reads/realigned.bam"
	rm "${OUTPUT}/bedfiles/${chr}/realigned_reads/paired_reads.bam" "${OUTPUT}/bedfiles/${chr}/realigned_reads/singletons.bam"

	#/usr/bin/time -v -p -o "${OUTPUT}/${chr}_airlift.bam.time" "${BINDIR}/samtools" merge --write-index -@ ${THREAD} --reference "${NEWREF}/${chr}.${SEQ_FILE_EXT}" -l5 -f "${OUTPUT}/${chr}_airlift.bam" "${OUTPUT}/bedfiles/${chr}/realigned_reads/paired_reads.bam" "${OUTPUT}/bedfiles/${chr}/realigned_reads/singletons.bam" "${OUTPUT}/bedfiles/${chr}/constant/constant_lifted.bam";
	/usr/bin/time -v -p -o "${OUTPUT}/${chr}_airlift.bam.time" "${BINDIR}/samtools" merge --write-index -@ ${THREAD} -l5 -f "${OUTPUT}/${chr}_airlift.bam" "${OUTPUT}/bedfiles/${chr}/constant/crossmap_final.bam" "${OUTPUT}/bedfiles/${chr}/realigned_reads/realigned.bam";

done

