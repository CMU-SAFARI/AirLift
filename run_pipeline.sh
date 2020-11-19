#!/bin/bash

SRCFOLDER=$1
#ref folders containing ref file:
OLDREF=$2
NEWREF=$3
SEQ_FILE_EXT=$4
READSIZE=$5
#make sure the bam files are indexed (i.e., samtools index)
READ_BAM=$6
FIRST_PAIR=$7
SECOND_PAIR=$8
OUTPUT=$9
THREAD=${10}

mkdir -p "${OUTPUT}/"

function generate_chain(){
    bash "${SRCFOLDER}/1-generate_chain/chain_install.sh" "${OUTPUT}/"
    bash "${SRCFOLDER}/1-generate_chain/sample_run.sh" \
    "${OLDREF}" \
    "${NEWREF}" \
    "${OUTPUT}/" \
    "${SRCFOLDER}/1-generate_chain/chain_generate.sh" \
    ${SEQ_FILE_EXT} \
    "${OUTPUT}/"
}

function generate_gaps(){
    for i in `echo ${OUTPUT}/*.chain`; do 
        local chr=`basename $i | sed s/.chain//`; \
        /usr/bin/time -v -p -o "${OUTPUT}/${chr}_extracted_gaps.fa.time" \
        python3 "${SRCFOLDER}/2-generate_gaps/extract_gaps.py" \
        $i \
        "${NEWREF}/$chr.${SEQ_FILE_EXT}" \
        ${READSIZE} \
        "${OUTPUT}/${chr}_extracted_gaps.fa"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`; \
        /usr/bin/time -v -p -o "${OUTPUT}/${chr}_kmer_gaps.fasta.time" \
        python3 "${SRCFOLDER}/2-generate_gaps/gaps_to_fasta.py" \
        "${OUTPUT}/${chr}_extracted_gaps.fa" \
        ${READSIZE} \
        "${OUTPUT}/${chr}_kmer_gaps.fasta" \
        10
    done
}

function align_gaps(){
    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`; \
        bash "${SRCFOLDER}/3-align_gaps/align_gaps.sh" \
        "${OLDREF}/${chr}.${SEQ_FILE_EXT}" \
        "${OUTPUT}/${chr}_kmer_gaps.fasta" \
        ${THREAD} \
        "${OUTPUT}/${chr}_aligned_gaps"
    done
}

function extract_reads(){
    mkdir "${OUTPUT}/bedfiles"
    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`; mkdir -p "${OUTPUT}/bedfiles/${chr}/"; \
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/extract_regions.time" \
        bash "${SRCFOLDER}/4-extract_reads/extract_regions.sh" \
        "${OUTPUT}/${chr}_aligned_gaps.bam" \
        "${OUTPUT}/bedfiles/${chr}/"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        local num_of_files=`ls -l ${OUTPUT}/bedfiles/${chr}/*.bed | wc -l`
        local sed_bed_files=$(echo ${OUTPUT}\/bedfiles\/${chr}/ | sed 's/\//\\\//g')
        local sed_read_scripts=$(echo ${SRCFOLDER}\/4-extract_reads\/extract_reads.bash | sed 's/\//\\\//g')
        local sed_read_bam=$(echo $READ_BAM | sed 's/\//\\\//g')
        cat "${SRCFOLDER}/4-extract_reads/slurm_extract_reads.sh" | \
        sed "s/{NUM_FILES}/$num_of_files/g" | \
        sed "s/{THREAD}/$THREAD/g" | \
        sed "s/{BED_FILES_FOLDER}/$sed_bed_files/g" | \
        sed "s/{EXTRACT_READS_SCRIPT}/$sed_read_scripts/g" | \
        sed "s/{READ_BAM_FILE}/$sed_read_bam/g" | \
        sed "s/{READSIZE}/$READSIZE/g" > "${OUTPUT}/${chr}_samplejob.sh"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        sbatch --wait "${OUTPUT}/${chr}_samplejob.sh"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed.time" \
        cat "${OUTPUT}/bedfiles/${chr}/"*.bed | sortBed -i - | mergeBed -i - > "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        mkdir -p "${OUTPUT}/bedfiles/${chr}/retired_bed/"
        mkdir -p "${OUTPUT}/bedfiles/${chr}/constant_bed/"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed.time" \
        python3 "${SRCFOLDER}/4-extract_reads/get_retired_regions.py" \
        ${READSIZE} \
        "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed" \
        $i \
        "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed.time" \
        bash "${SRCFOLDER}/4-extract_reads/extract_reads_noprune.sh" \
        ${READ_BAM} \
        "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed" > "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        mkdir -p "${OUTPUT}/bedfiles/${chr}/reads/"
        cat "${OUTPUT}/bedfiles/${chr}/"*.reads "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed" > "${OUTPUT}/bedfiles/${chr}/updated_and_retired_reads.bed"
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/extract_sequences_full.time" \
        bash "${SRCFOLDER}/4-extract_reads/extract_sequence.sh" \
        ${FIRST_PAIR} \
        ${SECOND_PAIR} \
        "${OUTPUT}/bedfiles/${chr}/updated_and_retired_reads.bed" \
        "${OUTPUT}/bedfiles/${chr}/reads/"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/reads/paired_reads.time" \
        bwa mem -M -t $THREAD "${NEWREF}/$chr.${SEQ_FILE_EXT}" "${OUTPUT}/bedfiles/${chr}/reads/reads_1.fastq" "${OUTPUT}/bedfiles/${chr}/reads/reads_2.fastq" | \
        samtools view -h -F4 | \
        samtools sort -m 16g -l0 > "${OUTPUT}/bedfiles/${chr}/reads/paired_reads.bam"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/reads/singletons_reads.time" \
        bwa mem -M -t $THREAD "${NEWREF}/$chr.${SEQ_FILE_EXT}" "${OUTPUT}/bedfiles/${chr}/reads/singletons.fastq" | \
        samtools view -h -F4 | \
        samtools sort -m 16g -l0 > "${OUTPUT}/bedfiles/${chr}/reads/singletons_reads.bam"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        bash "${SRCFOLDER}/0-align_reads.sh" \
        "${NEWREF}/$chr.${SEQ_FILE_EXT}" \
        "${OUTPUT}/bedfiles/${chr}/reads/reads" \
        "${OUTPUT}/bedfiles/${chr}/reads/paired" \
        $THREAD
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        bash "${SRCFOLDER}/0-align_singletons.sh" \
        "${NEWREF}/$chr.${SEQ_FILE_EXT}" \
        "${OUTPUT}/bedfiles/${chr}/reads/singletons.fastq" \
        "${OUTPUT}/bedfiles/${chr}/reads/singletons" \
        $THREAD
    done
}

function extract_reads_with_parallel(){
    mkdir "${OUTPUT}/bedfiles"
    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`; mkdir -p "${OUTPUT}/bedfiles/${chr}/"; \
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/extract_regions.time" \
        bash "${SRCFOLDER}/4-extract_reads/extract_regions.sh" \
        "${OUTPUT}/${chr}_aligned_gaps.bam" \
        "${OUTPUT}/bedfiles/${chr}/"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        local num_of_files=`ls -l ${OUTPUT}/bedfiles/${chr}/*.bed | wc -l`
        local sed_bed_files=$(echo ${OUTPUT}\/bedfiles\/${chr}/ | sed 's/\//\\\//g')
        local sed_read_scripts=$(echo ${SRCFOLDER}\/4-extract_reads\/extract_reads.bash | sed 's/\//\\\//g')
        local sed_read_bam=$(echo $READ_BAM | sed 's/\//\\\//g')
        cat "${SRCFOLDER}/4-extract_reads/parallel_extract_reads.sh" | \
        sed "s/{NUM_FILES}/$num_of_files/g" | \
        sed "s/{THREAD}/$THREAD/g" | \
        sed "s/{BED_FILES_FOLDER}/$sed_bed_files/g" | \
        sed "s/{EXTRACT_READS_SCRIPT}/$sed_read_scripts/g" | \
        sed "s/{READ_BAM_FILE}/$sed_read_bam/g" | \
        sed "s/{READSIZE}/$READSIZE/g" > "${OUTPUT}/${chr}_samplejob.sh"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        bash "${OUTPUT}/${chr}_samplejob.sh"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed.time" \
        cat "${OUTPUT}/bedfiles/${chr}/"*.bed | sortBed -i - | mergeBed -i - > "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        mkdir -p "${OUTPUT}/bedfiles/${chr}/retired_bed/"
        mkdir -p "${OUTPUT}/bedfiles/${chr}/constant_bed/"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed.time" \
        python3 "${SRCFOLDER}/4-extract_reads/get_retired_regions.py" \
        ${READSIZE} \
        "${OUTPUT}/bedfiles/${chr}/merged_${chr}.bed" \
        $i \
        "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed.time" \
        bash "${SRCFOLDER}/4-extract_reads/extract_reads_noprune.sh" \
        ${READ_BAM} \
        "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_regions.bed" > "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        mkdir -p "${OUTPUT}/bedfiles/${chr}/reads/"
        cat "${OUTPUT}/bedfiles/${chr}/"*.reads "${OUTPUT}/bedfiles/${chr}/retired_bed/retired_reads.bed" > "${OUTPUT}/bedfiles/${chr}/updated_and_retired_reads.bed"
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/extract_sequences_full.time" \
        bash "${SRCFOLDER}/4-extract_reads/extract_sequence.sh" \
        ${FIRST_PAIR} \
        ${SECOND_PAIR} \
        "${OUTPUT}/bedfiles/${chr}/updated_and_retired_reads.bed" \
        "${OUTPUT}/bedfiles/${chr}/reads/"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/reads/paired_reads.time" \
        bwa mem -M -t $THREAD "${NEWREF}/$chr.${SEQ_FILE_EXT}" "${OUTPUT}/bedfiles/${chr}/reads/reads_1.fastq" "${OUTPUT}/bedfiles/${chr}/reads/reads_2.fastq" | \
        samtools view -h -F4 | \
        samtools sort -m 16g -l0 > "${OUTPUT}/bedfiles/${chr}/reads/paired_reads.bam"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/reads/singletons_reads.time" \
        bwa mem -M -t $THREAD "${NEWREF}/$chr.${SEQ_FILE_EXT}" "${OUTPUT}/bedfiles/${chr}/reads/singletons.fastq" | \
        samtools view -h -F4 | \
        samtools sort -m 16g -l0 > "${OUTPUT}/bedfiles/${chr}/reads/singletons_reads.bam"
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        bash "${SRCFOLDER}/0-align_reads.sh" \
        "${NEWREF}/$chr.${SEQ_FILE_EXT}" \
        "${OUTPUT}/bedfiles/${chr}/reads/reads" \
        "${OUTPUT}/bedfiles/${chr}/reads/paired" \
        $THREAD
    done

    for i in `echo ${OUTPUT}/*.chain`; do
        local chr=`basename $i | sed s/.chain//`
        bash "${SRCFOLDER}/0-align_singletons.sh" \
        "${NEWREF}/$chr.${SEQ_FILE_EXT}" \
        "${OUTPUT}/bedfiles/${chr}/reads/singletons.fastq" \
        "${OUTPUT}/bedfiles/${chr}/reads/singletons" \
        $THREAD
    done
}


function main(){
    local chain_existed=$(ls ${OUTPUT}/*.chain 2> /dev/null)
    
    if [[ ${#chain_existed} -gt 0 ]]; then
        echo "Line "${LINENO}": In function "${FUNCNAME}": Use user offered chain "${chain_existed}
    else
        generate_chain
    fi

    generate_gaps

    align_gaps

    if [[ -z ${11} ]]; then
        extract_reads_with_parallel
    else
        extract_reads
    fi
}

#constant regions:
#for i in `echo ${OUTPUT}/*.chain`; do local chr=`basename $i | sed s/.chain//`; /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/constant_bed/constant_regions.bed.time" python3 "${SRCFOLDER}/4-extract_reads/get_constant_regions.py" ${READSIZE} $i "${OUTPUT}/bedfiles/${chr}/constant_bed/constant_regions.bed"; done
#sbatch --local wrap="for i in `echo ${OUTPUT}/*.chain`; do local chr=`basename $i | sed s/.chain//`; /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/constant_bed/constant_reads.bed.time" bash "${SRCFOLDER}/4-extract_reads/extract_reads_noprune.sh" ${READ_BAM} "${OUTPUT}/bedfiles/${chr}/constant_bed/constant_regions.bed" > "${OUTPUT}/bedfiles/${chr}/constant_bed/constant_reads.bed"; done"

