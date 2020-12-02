#!/bin/bash

BINDIR=$1
TARGET=$2 #old
QUERY=$3 #new
WORK_DIR=$4 #chain files will be created inside this directory

mkdir -p "${WORK_DIR}/tmp/"

"${BINDIR}/faSplit" -lift="${WORK_DIR}/tmp/query_sp.lft" size "${QUERY}" -oneFile 3000 "${WORK_DIR}/tmp/query_sp" #split new

"${BINDIR}/faToTwoBit" "${WORK_DIR}/tmp/query_sp.fa" "${WORK_DIR}/tmp/query_sp.2bit"
"${BINDIR}/twoBitInfo" "${WORK_DIR}/tmp/query_sp.2bit" stdout | sort -k2,2nr > "${WORK_DIR}/tmp/query_sp.sizes"

"${BINDIR}/faToTwoBit" "${QUERY}" "${WORK_DIR}/tmp/query.2bit" #generate 2bit without splitting.
"${BINDIR}/twoBitInfo" "${WORK_DIR}/tmp/query.2bit" stdout | sort -k2,2nr > "${WORK_DIR}/tmp/query.sizes"

"${BINDIR}/faToTwoBit" "${TARGET}" "${WORK_DIR}/tmp/target.2bit"
"${BINDIR}/twoBitInfo" "${WORK_DIR}/tmp/target.2bit" stdout | sort -k2,2nr > "${WORK_DIR}/tmp/target.sizes"

# twoBitToFa "${WORK_DIR}/tmp/target.2bit" stdout | faSize stdin #(the result is 2937639113)
# calc \( 2937639113 / 2861349177  \) \* 1024 #= 1050 (rounded down to the nearest 50 -> we use it below for repMatch)

mkdir "${WORK_DIR}/blat"
#generate ooc file based on the calculation above
"${BINDIR}/blat/blat" "${WORK_DIR}/tmp/target.2bit" /dev/null /dev/null -tileSize=11 -makeOoc="${WORK_DIR}/tmp/target.ooc"

"${BINDIR}/blat/blat" "${WORK_DIR}/tmp/target.2bit" "${WORK_DIR}/tmp/query_sp.2bit" "${WORK_DIR}/blat/target_to_query_sp.psl" -tileSize=11 -ooc="${WORK_DIR}/tmp/target.ooc" -minScore=100 -minIdentity=100 -fastMap

"${BINDIR}/liftUp" -pslQ "${WORK_DIR}/blat/target.psl" "${WORK_DIR}/tmp/query_sp.lft" warn "${WORK_DIR}/blat/target_to_query_sp.psl"

"${BINDIR}/axtChain" -linearGap=medium -psl "${WORK_DIR}/blat/target.psl" "${WORK_DIR}/tmp/target.2bit" "${WORK_DIR}/tmp/query.2bit" "${WORK_DIR}/target_to_query.chain"

rm -rf "${WORK_DIR}/tmp/"

