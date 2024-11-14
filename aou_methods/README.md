## Disclosure
The original analysis was done using VCFs from short read whole genome sequencing (v7) in All of Us (AoU) VCFs nor the associated covariate files can be exported from the AoU database. To provide interested parties with a ‘replica’ of our workflow we provide a snakemake pipeline which walks through the steps of our analysis. Follow the steps below to get started!

**The p-values, gene names, and variants output in this example are not correct. They are only meant to simulate our work on AoU. To ensure that there is adequate overlap between the available VCF from 1K Genomes and the example ClinVar data, the provided ClinVar table has been edited to include variants that are present in the VCF**

## Walk though of the pipeline
1)	Generate a genotype table from all VCF’s located in the designated input VCF directory. * A genotype table is a tab separated table with variant snp id’s in the first column followed by the genotype for each individual in subsequent columns *
2)	Assign inheritance to genes present in the table of ClinVar pathogenic/likely pathogenic variants, using OMIM. 
3)	For each disease/phenotype in your custom ClinVar table conduct logistic regression on the thyroid cancer diagnosis by sex, age, principal components 1-16, and whether that individual has at least one non reference allele for a known pathogenic or likely pathogenic variant associated with said disease in ClinVar. We report AIC and p-values for each model. 
4)	For each variant in your custom ClinVar table report the following:
a.	The associated disease from ClinVar
b.	The pathogenic variant ID (GRCh38)
c.	The total number of individuals with a positive thyroid cancer diagnosis and a pathogenic mutation
d.	A  pipe ‘|’ separated list of participant IID’s which were found to have both thyroid cancer and a mutation. 

## Additional scripts
In addition to what is required for the pipeline described above you will also find additional scripts used to confirm that the variants reported are in fact pathogenic or likely pathogenic and that the number of people with thyroid cancer and a mutation are correct. These were analyzed manually and cannot be successfully run outside AoU. 

Additional scripts:
* `aou_methods/R_scripts/APC_pathogenic_variants_AoU_v7.R`
* `aou_methods/R_scripts/PTEN_pathogenic_variants_AoU_v7.R`
* `aou_methods/R_scripts/RET_pathogenic_variants_AoU_v7.R`
* ‘`aou_methods/manually_check_pathogeniciy.py`

## Environment set-up 
```bash
git clone https://github.com/pozdeyevlab/clinvar-tc.git
cd clinvar-tc/aou_methods
conda env create -f environment.yml
conda activate aou_methods
poetry install
```
## Required downloads
•	This work requires that `bcftools` be installed on your machine. Please follow the instructions from [bcftools]( https://samtools.github.io/bcftools/) to install.
•	Please download the vcf for chromosome 21 from [1K Genomes]( https://www.internationalgenome.org/) by following the steps below: 
```bash
cd clinvar-tc/aou_methods
mkdir example_vcf
cd example_vcf
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz 
mv 1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz chr21.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
mv 1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi chr21.vcf.gz.tbi
```
•	OMIM Inheritance data can be downloaded by following the instructions on the OMIM website [here]( https://data.omim.org/downloads/). Make sure that the file `genemap2.txt` is  places in `example_inputs/`. 

## Launching the pipeline
```bash
cd aou_methods/
snakemake --cores 1 --configfile config.yaml
```

This will generate the following files in the default output directory
* dev_output/logistic_regression_and_overlap/variant_summary.tsv
* dev_output/logistic_regression_and_overlap/logistic_regression.tsv

`variant_summary` contains the following columns:
|Column Name    |Description     |
|---------------|----------------|
|disease |name of disease form ClinVar|
|variant |snp ID (GRCh38)|
|inheritance|inheritance from OMIM( genemap2.txt)|
|gene|gene name from ClinVar|
|sample_ids|pipe separated list of iids with both thyroid cancer and an alternate allele for the given variant|

`logistic_regression` contains the following columns:
|Column Name    |Description     |
|---------------|----------------|
|disease |name of disease form ClinVar|
|AIC |Akaike information criterion from the logistic regression model|
|p_value|P-value from the logistic regression model|
|total_non_ref_individuals|numbe rof individuals with at least one non reference allele for a variant associated with the given disease|
|cancer_and_disease|number of individuals with a mutation and thyroid cancer|
