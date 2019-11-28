sbatch -c 48 --wrap="./minimap2 -x sr -t 48 -a /mnt/panzer/firtinac/genomic_data/refs/hg38/hg38.fa /mnt/panzer/firtinac/genomic_data/na12878/illumina/platinum/ERR262997_1.fastq.gz /mnt/panzer/firtinac/genomic_data/na12878/illumina/platinum/ERR262997_2.fastq.gz | samtools sort -m 16G -l0 -o /mnt/panzer/firtinac/genomic_data/na12878/illumina/platinum/tmp.ERR262997.bam -"

