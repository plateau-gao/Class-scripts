nextflow run nf-core/rnaseq \
    --input GSE151090/sample.csv \
    --outdir s3://bl5632/students/plateau/results/GSE151090 \
    --genome GRCh37 \
    --awsqueue BL5632_Yuan_Zekun \
    --awsregion ap-southeast-1 \
    -profile docker,awsbatch \
    -work-dir s3://bl5632/students/plateau/work
