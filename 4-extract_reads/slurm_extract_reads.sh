#!/bin/bash
#
#SBATCH --job-name=ext_reads
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_ext_reads.out
#SBATCH --array=1-{NUM_FILES}%{THREAD}

#variables to change are:
#{BED_FILES_FOLDER}
#{READ_BAM_FILE}
#{EXTRACT_READS_SCRIPT}
#{NUM_FILES}
#{THREAD}

srun hostname
FILES=({BED_FILES_FOLDER}/*.bed)
BED_FILE=${FILES[${SLURM_ARRAY_TASK_ID}]}
OUT_READS=$(echo ${BED_FILE}.reads)

srun --output=${OUT_READS} \
/usr/bin/time -v -p -o "${OUT_READS}.time" \
bash {EXTRACT_READS_SCRIPT} \
{READ_BAM_FILE} \
${BED_FILE} \
{READSIZE}
sleep 1

