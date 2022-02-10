mkdir -p index
bowtie2-build data/sacCer3.fa,data/ty5_6p.fa index/sacCer3_mix
# align to reference
mkdir -p results/sam
bowtie2 -x index/sacCer3_mix -1 data/e0175724_1.fq -2 data/e0175724_2.fq -S results/sam/mix.sam &>results/sam/mix.txt
# filter reads
mkdir -p results/bam
samtools view -bS -F 14 results/sam/mix.sam > results/bam/discordant.bam
samtools view -bS -f 2 results/sam/mix.sam > results/bam/concordant.bam

samtools sort results/bam/discordant.bam -o results/bam/discordant-sort.bam
samtools sort results/bam/concordant.bam -o results/bam/concordant-sort.bam

samtools index results/bam/discordant-sort.bam
samtools index results/bam/concordant-sort.bam

bedtools genomecov -ibam results/bam/discordant-sort.bam -bg -max 10 | grep -v TY5 | awk '$4==10' | bedtools merge -d 1000 | awk '{printf "%s:%d-%d\n", $1, $2, $3}' > results/result.txt

