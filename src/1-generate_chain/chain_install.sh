#!/bin/bash

#requires an internet connection

UTIL_DIR=$1 #under this folder, this script will create a folder called "utils" in which the installation will happen.

mkdir -p "${UTIL_DIR}/utils"
chmod 755 "${UTIL_DIR}/utils"

case "$(uname -s)" in
	Darwin*) rsync_comm="rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/ ${UTIL_DIR}/utils/";;
	Linux*) rsync_comm="rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ${UTIL_DIR}/utils/";;
esac

rsync -aP $rsync_comm

curl -o "${UTIL_DIR}/utils/bedSingleCover.pl" 'http://genome-source.soe.ucsc.edu/gitlist/kent.git/raw/master/src/utils/bedSingleCover.pl'

