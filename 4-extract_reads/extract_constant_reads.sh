export SRCFOLDER=$1
export READSIZE=$2
export OLDREF=$3
export READ_BAM=$4
export NEWREF=$5
export OUTPUT=$6

for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; mkdir -p "${OUTPUT}/bedfiles/${chr}/constant_bed/"; done
sbatch --wrap="for i in `echo ${OUTPUT}/*.chain`; do chr=`basename $i | sed s/.chain//`; /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/constant_bed/constant_regions.bed.time" python "${SRCFOLDER}/4-extract_reads/get_constant_regions.py" ${READSIZE} $i "${OUTPUT}/bedfiles/${chr}/constant_bed/constant_regions.bed"; /usr/bin/time -v -p -o "${OUTPUT}/bedfiles/${chr}/constant_bed/constant_reads.bed.time" sh "${SRCFOLDER}/4-extract_reads/extract_reads.sh" ${READ_BAM} "${OUTPUT}/bedfiles/${chr}/constant_bed/constant_regions.bed" > "${OUTPUT}/bedfiles/${chr}/constant_bed/constant_reads.bed"; done"


