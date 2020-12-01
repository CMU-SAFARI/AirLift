# Extract the gaps from the genome FASTA file depending on the chain file and the size of the reads. 
python extract_masked_gaps.py <chain file> <genome FASTA file> <size of read> <output file>
#e.g., $ python extract_masked_gaps.py ../liftover_files/hg19ToHg38.over.chain ../hg38/hg38.fa <size of read> masked_output.txt

# Extract overlapping reads from the gap contigs extracted in script above. 
python gaps_to_fasta.py <input FASTA file> <read size> <output FASTA file> <len that each read skips by>


