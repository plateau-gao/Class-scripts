adapters=/usr/local/Caskroom/miniconda/base/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa
mkdir -p results/trimmed
mkdir -p results/orphanded
for file in data/*_1.fastq.gz
do
    SRR=$(basename $file _1.fastq.gz)
    echo working on $SRR
    trimmomatic PE data/${SRR}_1.fastq.gz data/${SRR}_2.fastq.gz \
        results/trimmed/${SRR}_1.trim.fastq.gz results/orphanded/${SRR}_1.untrim.fastq.gz \
        results/trimmed/${SRR}_2.trim.fastq.gz results/orphanded/${SRR}_2.untrim.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${adapters}:2:40:15
done
