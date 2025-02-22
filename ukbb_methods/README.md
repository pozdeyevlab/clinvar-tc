# UKBB Methods

## Data Used from UKBB
#### Genotypes were extracted from population exome VCF's released April 2022. These files are only available through the UKBB Research Analysis Platform (RAP). Exome data is included in UKBB bulk data release at cost tier 3.

## Code
#### For a working example please refer to the `aou_methods`. This code is meant to be used with dx-toolkit for use with UKBB. 

### Step 1: Querying Exome Files & Creating A Genotype Table
In total all 975 exome files were filtered for pathogenic variants listed in clinvar, and then compiled into one merged table with genotypes for each particiapnt. 

```bash
dx login
dx ls /path/to/exome/files/*.vcf.gz > file_list.tsv

# The following command will launch job for every file in `file_list.tsv`
bash bash_scripts/launch_bcftools.sh file_list.tsv

# Once all exomes have been queried create the header and concatenate the genotypes
# The easiest method is to launch a ttyd instance -- then follow these commands
dx download /path/to/any/exome
echo 'variant' > header.txt
bcftools query -l exome.vcf >> header.txt
sed -i.bak 's/\n/\t/g' header.txt
dx upload header.txt --path=/path/to/output/dir/

# Concat Genotypes
# Make sure that concat.sh from bash_scripts/concat.sh is available on UKBB in your input directory
bash concat_on_dx.sh

# Run regression
# In order to deal with the increased number of samples in UKBB, it is recommended to run each gene in clinvar seperately.

# Separate clinvar by gene, name, with gene and chromosome, and then upload to UKBB
# assumes that variant_summary_pathogenic_Dec_30_24.txt is in working directory
# assumes that a file genes.txt is in the working directoyr which contains a list of gnene symbols of interest
python simplify.py

bash launch_regression.sh genes.txt
```

This will output a variant level summary as well as a disease level association output. 





