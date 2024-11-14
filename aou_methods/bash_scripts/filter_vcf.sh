%%bash
# Filter VCFs for regions in the tsv
cd clinvar
head regions.tsv

input_dir="vcfs"
output_dir="filtered_vcfs"
regions_file="regions.tsv"
mkdir -p $output_dir

for input_vcf in "$input_dir"/*.vcf.bgz; do
    # Extract the file name without the path and extension
    vcf_base=$(basename "$input_vcf" .vcf.bgz)

    # Create the output file path
    output_vcf="$output_dir/${vcf_base}_filtered.vcf.gz"

    # Find number of available threads (multithreading only applicable for gzipping output vcfs)
    nproc=$(nproc)

    # Use bcftools to filter regions based on the regions.tsv file and split multi-allelic sites
    bcftools norm -m- -R "$regions_file" --threads "$nproc" -O z -o "$output_vcf" "$input_vcf"

    echo "Filtered $vcf_base to $output_vcf"
done

%%bash
# Filter VCFs for regions in the tsv
cd clinvar

input_dir="filtered_vcfs"
output_dir="genotype_tables"
mkdir -p $output_dir

for input_vcf in "$input_dir"/*.vcf.gz; do
    # Extract the file name without the path and extension
    vcf_base=$(basename "$input_vcf" .vcf.bgz)

    # Create the output file path
    output_tsv="$output_dir/${vcf_base}.tsv.gz"
    sample_tsv="$output_dir/${vcf_base}_sample.tsv.gz"

    # Find number of available threads (multithreading only applicable for gzipping output vcfs)
    nproc=$(nproc)

    # Per VCF make GT Table
    sample_ids=$(bcftools query -l $input_vcf | tr '\n' '\t')
    echo -e "variant\t$sample_ids" | bgzip -c >> $sample_tsv
    bcftools query -f 'chr%CHROM:%POS:%REF:%ALT[\t%GT]\n' $input_vcf| bgzip -c >>$output_tsv
    echo "GT Table Written"

done