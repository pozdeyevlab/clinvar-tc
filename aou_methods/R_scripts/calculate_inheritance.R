library(optparse)
# Check if the package is installed
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  # If not installed, install it from GitHub using remotes
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  # Install GenomicSEM from GitHub
  remotes::install_github("tidyverse/tidyverse")
}
library(data.table)
library(R.utils)
library(dplyr)
library(tidyverse)
#library(bigrquery)

# Define command-line arguments
option_list <- list(
  make_option(c("--omim"), type="character", help="genemap2.txt file which can be downloaded following instructions in the readme"),
  make_option(c("--clinvar"), type="character", help="Table containing pathogenic variants from clivar database"),
  make_option(c("--out"), type="character", help="Where to write the inheritance table")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Read in pathogenic clinvar table
clinvar <- fread(opt$clinvar)

# Split PhenotypeList into seperate rows on '|' (pipe) then on ';' if not surrounded in ()
clinvar_dt <- clinvar %>% separate_rows(PhenotypeList, sep = "\\|") %>%
    separate_rows(PhenotypeList, sep = ";(?![^()]*\\))", convert = TRUE)

# OMIM
# https://data.omim.org/downloads/kigz6AEgSeCNrRbGjgv1Hw/genemap2.txt
omim <- fread(opt$omim)

# Generate inheritance table
inheritance_df <- clinvar_dt %>% 
    dplyr::select(GeneSymbol, PhenotypeList) %>%
    separate_rows(GeneSymbol, sep = "\\;") %>%
    filter(GeneSymbol != 'not provided') %>%
    distinct()

inheritance_df <- merge(inheritance_df, omim[ ,c('Chromosome', 'Genomic Position Start', 'Genomic Position End', 'Approved Gene Symbol')], by.x='GeneSymbol', by.y='Approved Gene Symbol')

colnames(inheritance_df) <- c('GeneSymbol', 'PhenotypeList', 'Chromosome', 'Start', 'End')

inheritance_df <- inheritance_df[,1:5]
print(dim(inheritance_df))


for (index in 1:nrow(inheritance_df)) {
    gene <- inheritance_df$GeneSymbol[index]
    phenotype <- inheritance_df$PhenotypeList[index]
    
    from_name_recessive <- length(grep('autosomal recessive', tolower(phenotype)))
    from_name_dominant <- length(grep('autosomal dominant', tolower(phenotype)))
    
    if (!rlang::is_empty(gene)){
        sub_omim <- as.data.frame(omim[omim$'Approved Gene Symbol'==gene,] %>% 
                                  separate_rows(Phenotypes, sep = ";"))
        
        recessive <- length(grep('Autosomal recessive', sub_omim[grep(phenotype, sub_omim$Phenotypes),]))
        dominant <- length(grep('Autosomal dominant', sub_omim[grep(phenotype, sub_omim$Phenotypes),]))

        if (recessive==0 & dominant==1){
            inheritance <- 'dominant_from_omim'
        }
        if (recessive==1 & dominant==0){
            inheritance <- 'recessive_from_omim'
        }
        if (recessive==1 & dominant==1){
            inheritance <- 'found_both_treat_as_dominant'
        }
        if (recessive==0 & dominant==0){
            inheritance <- 'not_found_in_omim_treat_as_dominant'
        }
    }
    
    if (from_name_recessive==1){
        inheritance <- 'recessive_from_name'
    }
    if (from_name_dominant==1){
        inheritance <- 'dominant_from_name'
    }
    inheritance_df$inheritance[index] <- inheritance
}

write.table(inheritance_df, opt$out, sep='\t', quote=FALSE, row.names=FALSE)