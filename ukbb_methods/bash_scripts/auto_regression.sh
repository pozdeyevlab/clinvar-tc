#!/bin/sh

### Requirements: 
# dx-toolkit 
# dx login to RAP

### How to Run:
# Run this script using: 
# bash auto_regression.sh
# on your local command line

### Inputs:
data_dir="regression_feb_2025"
inheritance='inheritance.tsv'
clinvar=$1
covariates='covariates_filtered.tsv'
genotype='exome_with_zeros.tsv'
variantout=$2
lrout=$3
split="test_split.sh"
r_script='ukbb_regression_20250216.R'

### Goal:
# Run multiple_logistic_regresison.R on UKBB

### Define command to run on RAP (use as input to dx run with -icmd="{$regression}")
regression="bash ${split} ${r_script} ${inheritance} ${clinvar} ${covariates} ${genotype} ${variantout} ${lrout}"


### Submit to the RAP with dx run
dx run swiss-army-knife -iin="${data_dir}/${inheritance}" \
   -iin="${data_dir}/${clinvar}" \
   -iin="${data_dir}/${covariates}" \
   -iin="${data_dir}/${genotype}" \
   -iin="${data_dir}/${r_script}" \
   -iin="${data_dir}/${split}" \
   -icmd="${regression}" \
   --instance-type "mem2_ssd1_v2_x96" \
   --destination="regression_20250218_results" \
   --tag="${clinvar}" --name='regression' \
   --brief --yes --priority high \
   --ignore-reuse

