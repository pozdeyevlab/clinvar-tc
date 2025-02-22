input=$1
while IFS= read -r line
do
    tbi=$line.tbi
    
    txt=$(echo $line | sed "s/.vcf.gz/.txt/g")

    bash bcftools.sh $line $tbi $txt

done < $input