# All of Us Methods
## Discosure
The original analysis was done using VCFs from short read whole genome sequencing (v7) in All of Us (AoU) VCFs nor the associated covariate files can be exported from the AoU database. To provide interested parties with a ‘replica’ of our workflow we provide a snakemake pipeline which walks through the steps of our analysis. Follow the steps below to get started!

**The p-values, gene names, and variants output in this example are not correct. They are only meant to simulate our work on AoU**

## Walk though of the pipeline
1)	Generate a genotype table from all VCF’s located in the designated input VCF directory. * A genotype table is a tab separated table with variant snp id’s in the first column followed by the genotype for each individual in subsequent columns *
2)	Assign inheritance to genes present in the table of ClinVar pathogenic/likely pathogenic variants, using OMIM. 
3)	For each disease/phenotype in your custom ClinVar table conduct logistic regression on the thyroid cancer diagnosis by sex, age, principal components 1-16, and whether or not that individual has at least one non reference allele for a known pathogenic or likely pathogenic variant associated with said disease in ClinVar. We report AIC and p-values for each model. 
4)	For each variant in your custom ClinVar table report the following:
a.	The associated disease from ClinVar
b.	The pathogenic variant ID (GRCh38)
c.	The total number of individuals with a positive thyroid cancer diagnosis and a pathogenic mutation
d.	A  pipe ‘|’ separated list of participant IID’s which were found to have both thyroid cancer and a mutation. 

## Additional scripts
In addition to what is required for the pipeline described above you will also find additional scripts used to confirm that the variants reported are in fact pathogenic or likely pathogenic and that the number of people with thyroid cancer and a mutation are correct. These were analyzed manually and cannot be successfully run outside AoU. 

Additional scripts:
•	`aou_methods/R_scripts/APC_pathogenic_variants_AoU_v7.R`
•	`aou_methods/R_scripts/PTEN_pathogenic_variants_AoU_v7.R`
•	`aou_methods/R_scripts/RET_pathogenic_variants_AoU_v7.R`
•	‘`aou_methods/manually_check_pathogeniciy.py`

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
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.annotated.vcf.gz
mv 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.annotated.vcf.gz chr21.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.annotated.vcf.gz.tbi
mv 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.annotated.vcf.gz.tbi chr21.vcf.gz.tbi
