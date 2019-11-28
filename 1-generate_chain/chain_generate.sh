#!/bin/sh

TARGET=$1 #old
QUERY=$2 #new
UTIL_DIR=$3
WORK_DIR=$4 #chain files will be created inside this directory

export PATH="${UTIL_DIR}/utils":"${UTIL_DIR}/utils/blat":$PATH
mkdir -p "${WORK_DIR}/tmp/"

faSplit -lift="${WORK_DIR}/tmp/query_sp.lft" size "${QUERY}" -oneFile 3000 "${WORK_DIR}/tmp/query_sp" #split new

faToTwoBit "${WORK_DIR}/tmp/query_sp.fa" "${WORK_DIR}/tmp/query_sp.2bit"
twoBitInfo "${WORK_DIR}/tmp/query_sp.2bit" stdout | sort -k2,2nr > "${WORK_DIR}/tmp/query_sp.sizes"

faToTwoBit "${QUERY}" "${WORK_DIR}/tmp/query.2bit" #generate 2bit without splitting.
twoBitInfo "${WORK_DIR}/tmp/query.2bit" stdout | sort -k2,2nr > "${WORK_DIR}/tmp/query.sizes"

faToTwoBit "${TARGET}" "${WORK_DIR}/tmp/target.2bit"
twoBitInfo "${WORK_DIR}/tmp/target.2bit" stdout | sort -k2,2nr > "${WORK_DIR}/tmp/target.sizes"

# twoBitToFa "${WORK_DIR}/tmp/target.2bit" stdout | faSize stdin #(the result is 2937639113)
# calc \( 2937639113 / 2861349177  \) \* 1024 #= 1050 (rounded down to the nearest 50 -> we use it below for repMatch)

mkdir "${WORK_DIR}/blat"
#generate ooc file based on the calculation above
blat "${WORK_DIR}/tmp/target.2bit" /dev/null /dev/null -tileSize=11 -makeOoc="${WORK_DIR}/tmp/target.ooc"

blat "${WORK_DIR}/tmp/target.2bit" "${WORK_DIR}/tmp/query_sp.2bit" "${WORK_DIR}/blat/target_to_query_sp.psl" -tileSize=11 -ooc="${WORK_DIR}/tmp/target.ooc" -minScore=100 -minIdentity=100 -fastMap

liftUp -pslQ "${WORK_DIR}/blat/target.psl" "${WORK_DIR}/tmp/query_sp.lft" warn "${WORK_DIR}/blat/target_to_query_sp.psl"

axtChain -linearGap=medium -psl "${WORK_DIR}/blat/target.psl" "${WORK_DIR}/tmp/target.2bit" "${WORK_DIR}/tmp/query.2bit" "${WORK_DIR}/target_to_query.chain"

rm -rf "${WORK_DIR}/tmp/"

