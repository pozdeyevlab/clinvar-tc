#!/bin/sh

### Requirements: 
# dx-toolkit 
# dx login to RAP

### How to Run:
# Run this script using: 
# bash bcftools.sh
# on your local command line

### Inputs:
data_dir="regression_feb_2025"
inheritance='inheritance.tsv'
clinvar='small.txt'
covariates='covariates_and_diagnosis_20250215.tsv'
genotype='exact_subset_exome_20250216.tsv'
variantout='out_test_variant.tsv'
lrout='out_test_lr.tsv'
r_script='ukbb_regression_20250216.R'

### Goal:
# Run multiple_logistic_regresison.R on UKBB

### Define command to run on RAP (use as input to dx run with -icmd="{$regression}")
regression="Rscript ${r_script} --inheritance ${inheritance} --clinvar ${clinvar} --covariates ${covariates} --genotype ${genotype} --VariantOut ${variantout} --LROut ${lrout}"


### Submit to the RAP with dx run
dx run swiss-army-knife -iin="${data_dir}/${inheritance}" \
   -iin="${data_dir}/${clinvar}" \
   -iin="${data_dir}/${covariates}" \
   -iin="${data_dir}/${genotype}" \
   -iin="${data_dir}/${r_script}" \
   -icmd="${regression}" \
   --instance-type "mem2_ssd1_v2_x96" \
   --destination="regression_feb_2025" \
   --tag="regression" --name='regression' \
   --brief --yes --priority high \
   --ignore-reuse


### Expected run time & cost
# Time = ~40 minutes
# Cost = ~0.08 
