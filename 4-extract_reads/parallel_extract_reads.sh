#!/bin/bash
# Need GNU parallel

module load parallel
source $(which env_parallel.bash)

#variables to change are:
#{BED_FILES_FOLDER}
#{READ_BAM_FILE}
#{EXTRACT_READS_SCRIPT}
#{NUM_FILES}
#{THREAD}

FILES=({BED_FILES_FOLDER}/*.bed)
BED_FILE=${FILES[${SLURM_ARRAY_TASK_ID}]}
OUT_READS=$(echo ${BED_FILE}.reads)

# Check the parallel command first 
env_parallel -j3 -k --dryrun /usr/bin/time -v -p -o "${OUT_READS}.time" bash {EXTRACT_READS_SCRIPT} {READ_BAM_FILE} {1} {READSIZE} '>' {2} ::: $(echo ${FILES[@]}) ::: $(echo ${FILES[@]} | awk '{for(i=1;i<=NF;i++) printf $i".reads ";}')
env_parallel -j3 -k --progress /usr/bin/time -v -p -o "${OUT_READS}.time" bash {EXTRACT_READS_SCRIPT} {READ_BAM_FILE} {1} {READSIZE} '>' {2} ::: $(echo ${FILES[@]}) ::: $(echo ${FILES[@]} | awk '{for(i=1;i<=NF;i++) printf $i".reads ";}')

module unload parallel
