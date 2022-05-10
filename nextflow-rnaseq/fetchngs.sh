nextflow run nf-core/fetchngs \
    --input SRR_Acc_List.txt \
    --outdir s3://bl5632/students/plateau/results/test \
    --awsqueue BL5632_Yuan_Zekun \
    --awsregion ap-southeast-1 \
    -profile docker,awsbatch \
    -work-dir s3://bl5632/students/plateau/work
