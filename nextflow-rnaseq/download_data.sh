outdir="/Users/plateau/Desktop/5632_ca2/result/GSE151090/nf-core_result/salmon"
cd $outdir
for id in `seq 2 19`
do
    dir=`awk -F , -v num=$id 'NR==num {print $1}' ~/Desktop/5632_ca2/GSE151090/sample.csv | sed 's/\"//g'` 
    mkdir $dir
    aws s3 cp s3://bl5632/students/plateau/results/GSE151090/star_salmon/$dir/quant.sf $dir/quant.sf
    aws s3 cp s3://bl5632/students/plateau/results/GSE151090/star_salmon/$dir/quant.genes.sf $dir/quant.genes.sf
done