# first parameter is the bed file 
# second parameter is the fasta file which we extract sequences from 

while IFS=$'\t' read -ra tmp
do
  samtools faidx "$2" "${tmp[0]}:${tmp[1]}-${tmp[2]}"
done < "${1}"

