set -e
# stop the execution of a script if a command or pipeline has an error

SCRIPT_DIR=$(dirname $(realpath $0))
BASE_DIR=$(dirname $SCRIPT_DIR)

cd ${BASE_DIR}

echo working in ${pwd}

genoem=data/genome/ecoli_rel606.fasta

export BOWTIE2_INDEXES=${BASE_DIR}/data/genome

echo "building the index of Ec606"

bowtie2-build $genome data/genome/Ec606

mkdir -p results/sam results/bam results/bcf results/vcf 

for fq1 in results/trimmed/*_1.trim.fastq.gz
do
    SRR=$(basename $fq1 _1.trim.fastq.gz)
    echo "working with $SRR"
    fq1=results/trimmed/${SRR}_1.trim.fastq.gz
    fq2=results/trimmed/${SRR}_2.trim.fastq.gz
    sam=results/sam/${SRR}.sam
    bam=results/bam/${SRR}.bam
    sort_bam=results/bam/${SRR}_sorted.bam
    raw_bcf=results/bcf/${SRR}_raw.bcf
    variants=results/vcf/${SRR}_variants.vcf
    final_variants=results/vcf/${SRR}_final_variants.vcf
    bowtie2 -x Ec606 --very-fast -p 4 -1 ${fq1} -2 ${fq2} -S ${sam}
    samtools view -S -b ${sam} > ${bam}
    samtools sort ${bam} -o ${sort_bam}
    samtools index ${sort_bam}
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf
    vcfutils.pl varFilter $variants > final_variants
done

