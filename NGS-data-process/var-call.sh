mkdir -p results/bcf
for file in results/bam/*-sorted.bam
do
    SRR=$(basename $file -sorted.bam)
    echo running $SRR
    bcftools mpileup -O b -o results/bcf/${SRR}_raw.bcf \
        -f data/genome/ecoli_rel606.fasta results/bam/${SRR}-sorted.bam
    bcftools call --ploidy 1 -m -v -o results/bcf/${SRR}_variants.vcf results/bcf/${SRR}_raw.bcf
    vcfutils.pl varFilter results/bcf/${SRR}_variants.vcf > results/bcf/${SRR}_final_variants.vcf
done