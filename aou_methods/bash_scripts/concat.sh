#!/bin/bash

# Read in named command line args 
while getopts ":f:v:o:" opt; do
  case $opt in
    f) files="$OPTARG"
    ;;
    v) vcf="$OPTARG"
    ;;
    o) output="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1
    ;;
  esac
done

# Make header
echo $vcf
sample_ids=$(bcftools query -l $vcf | tr '\n' '\t')

echo -e "variant\t$sample_ids" | bgzip -c > $output

# Loop
echo $files
files=${files//,/' '}
for file in $files; 
do
    echo $file
    gunzip -c $file | bgzip -c >> $output
done
