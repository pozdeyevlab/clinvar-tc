input=$1
while read line;
do
    clin=$(ls *$line.txt)

    chrom=$(echo $clin | cut -d'_' -f1)

    varout=$line'_variantout.tsv'

    lrout=$line'_lrout.tsv'

    bash auto_regression.sh $clin $varout $lrout

done < $input