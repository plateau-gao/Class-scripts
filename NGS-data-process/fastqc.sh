mkdir -p results/fastqc
fastqc data/*.fastq* -o results/fastqc
cd results/fastqc
for filename in *.zip
do
    unzip $filename
done
cat */summary.txt > fastq_summaries.txt
