if (!require(optparse)) install.packages('optparse')
if (!require(data.table)) install.packages('data.table')
if (!require(R.utils)) install.packages('R.utils')
if (!require(dplyr)) install.packages('dplyr')
if (!require(rlang)) install.packages('rlang')
if (!require(tibble)) install.packages('tibble')
if (!require(tidyr)) install.packages('tidyr')
if (!require(bit64)) install.packages('bit64')

library(optparse)
library(data.table)
library(R.utils)
library(rlang)
library(tibble)
library(dplyr)
library(tidyr)
library(bit64)

# Define command-line arguments
option_list <- list(
  make_option(c("--inheritance"), type="character", help="File with inheritance data per disease & gene from OMIM"),
  make_option(c("--clinvar"), type="character", help="Table containing pathogenic variants from clivar database"),
  make_option(c("--covariates"), type="character", help="Covariantes containing sample_id, age, sex, race, PC's1-10 and thyroid cancer diagnosis"),
  make_option(c("--genotype"), type="character", help="Genotype table created in previous step"),
  make_option(c("--VariantOut"), type="character", help="Write variant level output"),
  make_option(c("--LROut"), type="character", help="Write logistic regression output")
)
opt <- parse_args(OptionParser(option_list=option_list))



# Read in clinvar & inheritance data & covariates & genotypes
clinvar_dt <- fread(opt$clinvar)
clinvar_dt <- data.frame(clinvar_dt) %>% mutate(ref_first_variant = paste0('chr',Chromosome,':',PositionVCF,':',ReferenceAlleleVCF,':',AlternateAlleleVCF))
inheritance_df <- fread(opt$inheritance, sep='\t')

clinvar_dt <- clinvar_dt %>% separate_rows(PhenotypeList, sep = "\\|") %>%
    separate_rows(PhenotypeList, sep = ";(?![^()]*\\))", convert = TRUE)
clinvar_df = clinvar_dt
write.table('test.tsv', opt$out, sep='\t', quote=FALSE, row.names=FALSE)

covar <- fread(opt$covariates)
covar$sex_at_birth <- as.numeric(as.factor(covar$sex_at_birth))
#covar$race <- as.numeric(as.factor(covar$race))

# Read in genotype data
gt_table <- fread(opt$genotype, sep='\t', header=TRUE)

# Handle known fread error
gt_table <- gt_table[,1:(ncol(gt_table)-1)]


# For the purposes of this example this is a made-up list, in our actual analysis this is a list of meddullary thyroid cancer cases
samples_to_remove <- c('HG00096')


### OUTPUT CLASSES ###
# Set classes to record results

setClass(Class = "GLM",
  representation(
    AIC = "numeric",
    pvalue = "numeric",
    total = "numeric"
  )
)

setClass(Class = "Inheritance",
  representation(
    df = "data.frame"
  )
)

setClass(Class = "Inputs",
  representation(
    covariants = "data.frame",
    clinvar = "data.frame",
    inheritance = "data.frame",
    mt_samples = "list"
  )
)


setClass(Class = "MainResults",
  representation(
    AIC = "numeric",
    pvalue = "numeric",
    overlap_count = "numeric",
    total = "numeric"
  )
)

setClass(Class = "FinalTables",
  representation(
    variant_df = "data.frame",
    covariate_df = "data.frame"))


### WORKER FUNCTIONS ###
generate_output_table <- function(clinvar_dt) {
  output_df <- data.frame("disease" = unique(clinvar_dt$PhenotypeList))
  # Add columns for adj_r_squared, and p_value
  output_df$AIC <- rep(".", nrow(output_df))
  output_df$p_value <- rep(".", nrow(output_df))
  output_df$total_non_ref_individuals <- rep(".", nrow(output_df))
  output_df$cancer_and_disease <- rep(".", nrow(output_df))
  return(output_df)
}


generate_variant_output_table <- function() {
  variant_df <- data.frame(
    phenotype = character(0),
    variant = character(0),
    inheritance = character(0),
    sample_ids = character(0)
  )
  return(variant_df)
}


# Generate Variant List
generate_varlist <- function(clinvar_table, phenotype) {
  gene_symbol_df <- clinvar_table[clinvar_table$PhenotypeList == phenotype, ]

  varlist <- paste0(
    "chr",
    gene_symbol_df$Chromosome,
    ":",
    as.character(gene_symbol_df$PositionVCF),
    ":",
    as.character(gene_symbol_df$ReferenceAlleleVCF),
    ":",
    as.character(gene_symbol_df$AlternateAlleleVCF)
  )
  return(varlist)
}


# Filter genotype table according to varlist
filter_gt_dt <- function(gt_dt, varlist) {
  return(gt_dt[gt_dt$variant %in% varlist, ])
}


# Find inheritance according to variant, gene, and phenotype
find_inheritance <- function(inheritance_df, phenotype, clinvar_df, variant) {
  # Search inheritance table for inheritance
  # everything except diseases found to be explicitly recessive will be treated as dominant
  genes <- as.list(unique(clinvar_df[clinvar_df$ref_first_variant == variant, ]$GeneSymbol))
  for (gene in genes) {
    if (!rlang::is_empty(gene)) {
      if (!grepl("none of which curated to show dosage sensitivity", gene)) {
        inheritance <- inheritance_df[inheritance_df$GeneSymbol == gene &   inheritance_df$PhenotypeList == phenotype, ]$inheritance
        if (!rlang::is_empty(inheritance)) {
          return(inheritance)
        } else {
          return("not_found_in_omim_treat_as_dominant")
        }
      }
    } else {
      return("not_found_in_omim_treat_as_dominant")
    }
  }
}


# Calculate non reference variants according to inheritance
apply_inheritance_logic <- function(filtered_gt_dt, phenotype, inputs) {
  master <- as.data.frame(filtered_gt_dt[, 'variant'])
  # Empty column for inheritance
  master$inheritance <- rep("", nrow(master))
  for (index in 1:nrow(master)) { 
    variant <- master$variant[index]
    inheritance <- find_inheritance(inputs@inheritance,phenotype,inputs@clinvar,variant) 
    master$inheritance[index] <- inheritance
  }
  #
  dominant_rows <- which(grepl("dominant", master$inheritance))
  recessive_rows <- which(grepl("recessive", master$inheritance))

  if (!rlang::is_empty(dominant_rows) && !rlang::is_empty(recessive_rows)) {
    if (length(dominant_rows) == 1){
        dominant_df <- as.data.frame(t(apply(filtered_gt_dt[dominant_rows, -c('variant')], 2,
                          function(col) ifelse(col == "1/0" | col == "1/1" | col == '0/1', 1, 0))))
    } else {
        dominant_df <- as.data.frame(apply(filtered_gt_dt[dominant_rows, -c('variant')], 2,
                          function(col) ifelse(col == "1/0" | col == "1/1" | col == '0/1', 1, 0)))
    } 
    if (length(recessive_rows) == 1) {
        recessive_df <- as.data.frame(t(apply(filtered_gt_dt[recessive_rows, -c('variant')], 2,
                          function(col) ifelse(col == "1/1", 1, 0))))
    }
    else {
        recessive_df <- as.data.frame(apply(filtered_gt_dt[recessive_rows, -c('variant')], 2,
                          function(col) ifelse(col == "1/1", 1, 0)))  
    }
    dominant_df$variant_id <- filtered_gt_dt[dominant_rows, variant]                         
    recessive_df$variant_id <- filtered_gt_dt[recessive_rows, variant]
    dominant_df$inheritance <- master[dominant_rows, ]$inheritance
    recessive_df$inheritance <- master[recessive_rows, ]$inheritance
    final <- rbind(dominant_df, recessive_df)
    rm(dominant_df, recessive_df)                       
    return(final)
  }
  if (!is_empty(dominant_rows) && is_empty(recessive_rows)){
    if (length(dominant_rows) == 1){
        dominant_df <- as.data.frame(t(apply(filtered_gt_dt[dominant_rows, -c('variant')], 2,
                          function(col) ifelse(col == "1/0" | col == "1/1" | col == '0/1', 1, 0))))
    } else {
        dominant_df <- as.data.frame(apply(filtered_gt_dt[dominant_rows, -c('variant')], 2,
                          function(col) ifelse(col == "1/0" | col == "1/1" | col == '0/1', 1, 0)))
    } 
    dominant_df$variant_id <- filtered_gt_dt[dominant_rows, variant]
    dominant_df$inheritance <- master[dominant_rows, ]$inheritance
    return(dominant_df)
  }
  if (is_empty(dominant_rows) && !is_empty(recessive_rows)){
    if (length(recessive_rows) == 1) {
        recessive_df <- as.data.frame(t(apply(filtered_gt_dt[recessive_rows, -c('variant')], 2,
                          function(col) ifelse(col == "1/1", 1, 0))))
    }
    else {
        recessive_df <- as.data.frame(apply(filtered_gt_dt[recessive_rows, -c('variant')], 2,
                          function(col) ifelse(col == "1/1", 1, 0)))  
    }
    recessive_df$variant_id <- filtered_gt_dt[recessive_rows, variant]
    recessive_df$inheritance <- master[recessive_rows, ]$inheritance
    return(recessive_df)
  }
}

                       
write_secondary_output <- function(genotype_df, phenotype, clinvar_df){
    non_zero_index <- which(colSums(genotype_df[,!names(genotype_df) %in% c("variant_id", "inheritance")])>0)
    if (length(non_zero_index) == 1){
        colname <- names(non_zero_index)
        gt_only <- as.data.frame(genotype_df[,non_zero_index])
        colnames(gt_only) <- colname
        
    } else {
        gt_only <- as.data.frame(genotype_df[,non_zero_index])
    }
    for (index in 1:nrow(genotype_df)){
        inheritance <- genotype_df$inheritance[index]
        variant <- genotype_df$variant_id[index]
        if (sum(gt_only[index, ]) > 0) {
            tmp_index_df <- t(gt_only[index, ])
            non_zero <- which(rowSums(tmp_index_df)>0)
            tmp_index_df <- rownames_to_column(as.data.frame(tmp_index_df))
            samples <- (tmp_index_df[non_zero,]$rowname)
            collapsed_samples <- paste0(samples, collapse = "|")
            gene <- clinvar_df[clinvar_df$PhenotypeList == phenotype & clinvar_df$ref_first_variant == variant,]$GeneSymbol
            print(gene)
            tmp <- data.frame(disease=phenotype,
                               variant=variant,
                               inheritance=inheritance,
                               genename=gene,
                               cancer_and_mutation_ids=collapsed_samples)
            print(tmp)
            write.table(tmp, file=opt$VariantOut, append=TRUE, sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
        }
    }
}
                                       
                       
# Generate table for logistic regression
generate_logistic_model_df <- function(table, inputs, phenotype) {
  covar <- inputs@covariants
  men <- grepl('Multiple endocrine neoplasia', phenotype)
  mt <- grepl('Medullary thyroid ', phenotype)
  if (!men & !mt){
    # Remove medullary cases from model's not looking at Medullary Thyroid Cancer and MEN
    table <- table[,!names(table) %in% inputs@mt_samples]
    cols_to_keep <- intersect(names(table), covar$person_id)
    cols_to_keep_2 <- append(cols_to_keep, 'inheritance')
    cols_to_keep_3 <- append(cols_to_keep_2, 'variant_id')
    table <- table[, cols_to_keep_3, drop = FALSE]
    covar <- covar[covar$person_id %in% cols_to_keep_3,]
    covar$variant_count <- colSums(table[,!names(table) %in% c("variant_id", "inheritance")])
    tc_samples <- covar[covar$tc_diagnosis == 1, ]$person_id
    covar$variant_count <- colSums(
      table[ ,!names(table) %in% c("variant_id", "inheritance")]) 
    covar <- covar %>%
      mutate(variant_count = ifelse(variant_count >= 1, 1, 0)) %>%
      dplyr::select(-person_id)
    
    final <- table[,names(table) %in% c("variant_id", "inheritance", tc_samples)]
    return(new("FinalTables",
               variant_df=final,
               covariate_df=covar))
  } else {
    cols_to_keep <- intersect(names(table), covar$person_id)
    cols_to_keep_2 <- append(cols_to_keep, 'inheritance')
    cols_to_keep_3 <- append(cols_to_keep_2, 'variant_id')
    table <- table[, cols_to_keep_3, drop = FALSE]
    covar <- covar[covar$person_id %in% cols_to_keep_3,]
    covar$variant_count <- colSums(table[,!names(table) %in% c("variant_id", "inheritance")])
    tc_samples <- covar[covar$tc_diagnosis == 1, ]$person_id
    covar <- covar %>%
      mutate(variant_count = ifelse(variant_count >= 1, 1, 0)) %>%
      dplyr::select(-person_id)
    final <- table[,names(table) %in% c("variant_id", "inheritance", tc_samples)]
    return(new("FinalTables",
               variant_df=final,
               covariate_df=covar))
  }
}

                       
# Compelte logistic regression
conduct_logistic_regression <- function(log_df) {
    # Assign variant as binary column 0/1 for presence of atleast one disease associated variant
    # Calculate overlap betweeen cancer and disease associated variant
    # Calculate total number of individuals with disease associated variants
    model <- glm(tc_diagnosis ~ ., data = log_df, family = binomial)
    aic <- AIC(model)
    summary <- summary(model)
    p_value <- summary$coefficients['variant_count', "Pr(>|z|)"]
    total_variant <- sum(log_df$variant_count)
    return(new("GLM",
            AIC=aic,
            pvalue=p_value,
            total=total_variant))
}


# Main funciton
main <- function(inputs, gt_dt, phenotype, clinvar_df) {
    # Subset clinvar table and create varlist
    varlist <- generate_varlist(inputs@clinvar, phenotype)

    filtered_gt_dt <- gt_dt[gt_dt$variant %in% as.list(varlist),]
    # If there are 0 rows in the filtered table, then associated variants are not present in AoU dataset v7
    rows <- nrow(filtered_gt_dt)
    if (rows >= 1) {
        # For each row in filtered_gt_dt calculate disease presence based on available inheritance data
        # Table of IID's and the whether they are 0 (ref) or (1) non-ref at each variant (accounting for heritability)
        inheritance_results <- apply_inheritance_logic(filtered_gt_dt, phenotype, inputs)
        print('Successfully collected inheritance')
        
        # Format table so that it is compatible with downstream analysis
        # Generate data frame for linear regression
        final_tables <- generate_logistic_model_df(inheritance_results, inputs, phenotype)

        ## Calculate TC mutation overlap 
        cancer_and_disease_count <- nrow(final_tables@covariate_df %>% filter(tc_diagnosis == 1 & variant_count >= 1))
	print(cancer_and_disease_count)
        if (cancer_and_disease_count >= 1){
            # Conduct regression
            glm_results <- conduct_logistic_regression(final_tables@covariate_df)
            print('Successfully made logistic regression model')
            
            # Write details to second file
            write_secondary_output(final_tables@variant_df, phenotype, clinvar_df)
            
            ## Make and return instance of MainResults
            return(new("MainResults",
                AIC=glm_results@AIC,
                pvalue=glm_results@pvalue,
                overlap_count=cancer_and_disease_count,
                total=glm_results@total))
        } else {
            return(new("MainResults",
                AIC=NaN,
                pvalue=NaN,
                overlap_count=cancer_and_disease_count,
                total=sum(final_tables@covariate_df$variant_count)))
        }
    } else {
        return(new("MainResults",
            AIC=NaN,
            pvalue=NaN,
            overlap_count=NaN,
            total=NaN))
    }
}

### RUN WORK HERE ###
output_table <- generate_output_table(clinvar_dt)

variant_level_df <- generate_variant_output_table()
write.table(variant_level_df, file=opt$LROut, sep='\t', row.names=FALSE, quote=FALSE)

inputs <- new("Inputs",
    covariants=covar,
    clinvar=clinvar_df,
    inheritance=inheritance_df,
    mt_samples=as.list(samples_to_remove[1])
   )

# Subset to the first 3 phenotypes for demonstration purposes
#output_table = output_table[1:5,]

### CALL MAIN COMMAND ###
print(Sys.time())
for (index in 1:nrow(output_table)) {
    # Collect phenptype
    phenotype <- output_table$disease[index]
    # Call main function
    log_results <- main(inputs, gt_table, phenotype, clinvar_df) 
    output_table$AIC[index] <- log_results@AIC
    output_table$p_value[index] <- log_results@pvalue
    output_table$total_non_ref_individuals[index] <- log_results@total
    output_table$cancer_and_disease[index] <- log_results@overlap_count
    write.table(output_table, file=opt$LROut, sep='\t', row.names=FALSE, quote=FALSE)
}

# Add header to variant level summary
if (file.exists(opt$VariantOut)) {
  variant_results <- read.table(opt$VariantOut, sep='\t')
  variant_results <- unique(variant_results)
  colnames(variant_results) <- c('disease', 'variant','inheritance', 'gene', 'sample_ids')
  write.table(variant_results, file = opt$VariantOut, sep='\t', quote=FALSE, row.name=FALSE)
}
print(Sys.time())
