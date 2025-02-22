#!/bin/sh

### Requirements: 
# dx-toolkit 
# dx login to RAP

### How to Run:
# Run this script using: 
# bash bcftools.sh
# on your local command line

### Inputs:
data_input_vcf=$1
data_input_vcf_index=$2
data_output=$3

### Goal:
# Use bcftools to query the exomes for specific regions

# input and output directories
data_dir="Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release"
clinvar="clinvar_regions.tsv"

### Define command to run on RAP (use as input to dx run with -icmd="{$run_bcftools}")
run_bcftools="bcftools norm -m- -R ${clinvar} ${data_input_vcf} |\
 bcftools query -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n' > ${data_output}"

### Submit to the RAP with dx run
dx run swiss-army-knife -iin="${data_dir}/${data_input_vcf}" \
   -iin="${data_dir}/${data_input_vcf_index}" \
   -iin="SW_TEST/${clinvar}" \
   -icmd="${run_bcftools}" \
   --instance-type "mem1_ssd1_v2_x4" \
   --destination="SW_TEST" \
   --tag="Bcftools" --name='bcftools' \
   --brief --yes --priority high \
   --ignore-reuse


### Expected run time & cost
# Time = ~40 minutes
# Cost = ~0.08 
