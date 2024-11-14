install.packages('R.utils')
install.packages("eeptools")
library(data.table)
library(eeptools)

all <- fread("./PHENOTYPES/srWGS_all_demographics.csv")
tc <- fread("./PHENOTYPES/thyroid_cancer_srWGS.csv") 
clinical <- fread("./PHENOTYPES/srWGS_with_clinical_data_EHR_or_survey_demographics.csv")
pheo <- fread("./PHENOTYPES/pheochromocytoma_srWGS.csv")
phpt <- fread("./PHENOTYPES/primary_hyperparathyroidism_srWGS.csv")
age_at_diag <- fread("./PHENOTYPES/thyroid_cancer_age_at_diagnosis.csv")

ret = as.matrix(fread("./GENOTYPES/RET_genotypes.tsv.gz"))
dim(ret)

# last column has NAs removing
ret <- ret[, 1:(ncol(ret)-1)]

# getting the list of pathogenic variants from ClinVar table
vars <- fread("./GENES/RET.tsv.gz")

# pathalleles <- c(paste0("chr", 10, ":", vars$Start, ":", vars$ReferenceAlleleVCF, ":", vars$AlternateAlleleVCF),
#                 paste0("chr", 10, ":", vars$Start, ":", vars$AlternateAlleleVCF, ":", vars$ReferenceAlleleVCF))
# 
# pathalleles <- c(paste0("chr", 10, ":", vars$Start, ":", vars$ReferenceAlleleVCF, ":", vars$AlternateAlleleVCF))
#                  
# 
# # manually curated variants against the pull for RET gene from clinVar
# # no issues with strand flips
# ret[, "variant"][!(ret[, "variant"] %in% pathalleles)]
# ret[, "variant"][ret[, "variant"] %in% pathalleles]




# manually curated for variants that have the same locus but different alternate allele (to account for novel amino acid substitutions in 
# known hot spots)
keep <- c("chr10:43114598:G:C", "chr10:43114598:G:T", "chr10:43119548:G:A", "chr10:43100614:C:T", "chr10:43113622:G:A", "chr10:43113622:G:T", 
          "chr10:43113628:G:A", "chr10:43113648:T:A", "chr10:43113648:T:C", "chr10:43113649:G:C", "chr10:43113654:T:A", "chr10:43113654:T:G", 
          "chr10:43113655:G:C", "chr10:43114500:T:C", "chr10:43114500:T:G", "chr10:43114501:G:A", "chr10:43118392:G:C", "chr10:43119548:G:C", 
          "chr10:43119548:G:T", "chr10:43120144:T:G", "chr10:43120184:C:T")

# filtering for manually curated P/LP variants
ret <- ret[ret[, "variant"] %in% keep, ]
dim(ret)

# removing non-genotype column
varnames <- ret[, "variant"]

ret <- ret[,!colnames(ret) %in% c("variant")]
rownames(ret) <- varnames

# identifying individuals with P/LP variants
muts <- data.frame(variants = varnames, 
                   counts = numeric(length(varnames)))

mutants <- ret[, ret[1,] != "0/0" & ret[1,] != "./.", drop = F]
muts$counts[1] <- dim(mutants)[2]

for (i in 2:nrow(ret)){
  
  tmp <- ret[, ret[i,] != "0/0" & ret[i,] != "./.", drop = F]  
  print(dim(tmp))
  mutants <- cbind(mutants, tmp)
  
  muts$counts[i] <- dim(tmp)[2]
  
}


# adding annotations from MTC guidelines (Wells 2015) and ClinVar

muts

muts$manualannotation <- character(nrow(ret))

muts$manualannotation[1] <- "R77C, one star, pathogenic, MEN2, makes full length protein but affects glycosilation, not in guidelines"
muts$manualannotation[2] <- "C609W, one star, pathogenic, condition not provided, functional studies indicate pathogenicity, other amino acid substitution in this position are in guidelines"
muts$manualannotation[3] <- "C609F, two stars, pathogenic, MEN2, IN GUIDELINES, moderate risk"
muts$manualannotation[4] <- "C611Y, two stars, pathogenic, MEN2, IN GUIDELINES, moderate risk"
muts$manualannotation[5] <- "C618S, two starts, pathogenic, MEN2, IN GUIDLINES, moderate risk"
muts$manualannotation[6] <- "C6186, two starts, pathogenic, MEN2, IN GUIDLINES, moderate risk"
muts$manualannotation[7] <- "C618S, two starts, pathogenic, MEN2, IN GUIDLINES, moderate risk"
muts$manualannotation[8] <- "C620S, two starts, pathogenic, MEN2, IN GUIDLINES, moderate risk"
muts$manualannotation[9] <- "C620G, two starts, pathogenic, MEN2, not in the guidelines"
muts$manualannotation[10] <- "C620S, two starts, pathogenic, MEN2, IN GUIDLINES, moderate risk"
muts$manualannotation[11] <- "C634R, two starts, pathogenic, MEN2, IN GUIDELINES, high risk"
muts$manualannotation[12] <- "C634G, two starts, pathogenic, MEN2, IN GUIDELINES, high risk"
muts$manualannotation[13] <- "C634Y, two starts, pathogenic, MEN2, IN GUIDELINES, high risk"
muts$manualannotation[14] <- "K666N, two starts, pathogenic, MEN2, not in the guidelines but variant in the same amino acid is in the guidelines"
muts$manualannotation[15] <- "K666N, two starts, pathogenic, MEN2, not in the guidelines but variant in the same amino acid is in the guidelines"
muts$manualannotation[16] <- "E768D, two starts, pathogenic, MEN2, IN GUIDELINES, moderate risk"
muts$manualannotation[17] <- "V804M, two starts, pathogenic, MEN2, IN GUIDELINES, moderate risk"
muts$manualannotation[18] <- "V804L, two starts, pathogenic, MEN2, IN GUIDELINES, moderate risk"
muts$manualannotation[19] <- "V804L, two starts, pathogenic, MEN2, IN GUIDELINES, moderate risk"
muts$manualannotation[20] <- "S891A, two starts, pathogenic, MEN2, IN GUIDELINES, moderate risk"
muts$manualannotation[21] <- "S904F, two starts, likely pathogenic, not in the guidelines"


# 91 pathogenic variants that lead to mutations reported in guidelines. (counted based on the manual annotation)

write.csv(mutants, "./GENOTYPES/patients_with_pathogenic_mutations_in_RET.csv")
write.csv(muts, "./GENOTYPES/COUNTS_of_pathogenic_mutations_in_RET.csv")

# filtering

as.numeric(colnames(mutants)) %in% tc$person_id
as.numeric(colnames(mutants)) %in% phpt$pheo
as.numeric(colnames(mutants)) %in% phpt$person_id


# demographics for all APC mutants
all_RET <- all[all$person_id %in% as.numeric(colnames(mutants)), ]
ages <- age_calc(as.Date(all_RET$date_of_birth), enddate = Sys.Date(), units = "years")
mean(ages)
sd(ages)
table(all_RET$sex_at_birth)

# looking at individuals with thyroid cancer diagnosis and calculating demographics
ids <- as.numeric(colnames(mutants)[as.numeric(colnames(mutants)) %in% tc$person_id])
aad <- age_at_diag[age_at_diag$person_id %in% ids]
mean(aad$age_at_diagnosis)
sd(aad$age_at_diagnosis)
