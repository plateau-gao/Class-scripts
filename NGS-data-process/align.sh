mkdir -p results/sam
for file in results/trimmed/*1.trim.fastq.gz
do
    SRR=$(basename $file _1.trim.fastq.gz)
    echo running $SRR
    bowtie2 --very-fast -p 4 \
        -x data/genome/ec606 \
        -1 results/trimmed/${SRR}_1.trim.fastq.gz -2 results/trimmed/${SRR}_2.trim.fastq.gz \
        -S results/sam/${SRR}.sam
done