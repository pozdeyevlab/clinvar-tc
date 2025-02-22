#!/bin/sh

### Requirements: 
# dx-toolkit 
# dx login to RAP

### How to Run:
# Run this script using: 
# bash concat.sh
# on your local command line

### Goal:
# Concatenate all clinvar variants

# input and output directories
data_dir="TEST"
out_dir="working_dir_syndromic"
file_list="in_dir.tsv"
concat_helper_script="concat.sh"
output="concat_clinvar_exomes_jan_8_2025.tsv"

### Define command to run on RAP (use as input to dx run with -icmd="{$concat}")
concat="cp /mnt/project/${data_dir}/ukb23157_* . ;\
   bash ${concat_helper_script} -f ${file_list} -o ${output}"

### Submit to the RAP with dx run
dx run swiss-army-knife -iin="${data_dir}/${file_list}" \
   -iin="${data_dir}/${concat_helper_script}" \
   -icmd="${concat}" \
   --instance-type "mem1_ssd2_v2_x8" \
   --destination="working_dir_syndromic" \
   --tag="concat_test" --name='concat_test' \
   --brief --yes --priority high \
   --ignore-reuse


### Expected run time & cost
# Time = 
# Cost = 
