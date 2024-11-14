install.packages('R.utils')
install.packages("eeptools")
library(data.table)
library(eeptools)


ret = as.matrix(fread("./GENOTYPES/APC_genotypes.tsv.gz"))
dim(ret)

all <- fread("./PHENOTYPES/srWGS_all_demographics.csv")
tc <- fread("./PHENOTYPES/thyroid_cancer_srWGS.csv") 

age_at_diag <- fread("./PHENOTYPES/thyroid_cancer_age_at_diagnosis.csv")

# last column has NAs removing

ret <- ret[, 1:(ncol(ret)-1)]

# getting the list of pathogenic variants from ClinVar table
vars <- fread("./GENES/APC.tsv.gz")

pathalleles <- c(paste0("chr", vars$Chromosome[1], ":", vars$Start, ":", vars$ReferenceAlleleVCF, ":", vars$AlternateAlleleVCF))


# ret[, "variant"][!(ret[, "variant"] %in% pathalleles)]
keep <- ret[, "variant"][ret[, "variant"] %in% pathalleles]

 
# filtering for P/LP variants that matched
ret <- ret[ret[, "variant"] %in% keep, ]

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
  mutants <- cbind(mutants, tmp)
  muts$counts[i] <- dim(tmp)[2]
}

print(muts)
print(mutants)


write.csv(mutants, "./GENOTYPES/patients_with_pathogenic_mutations_in_APC.csv")
write.csv(muts, "./GENOTYPES/COUNTS_of_pathogenic_mutations_in_APC.csv")

# demographics for all APC mutants
all_APC <- all[all$person_id %in% as.numeric(colnames(mutants)), ]
ages <- age_calc(as.Date(all_APC$date_of_birth), enddate = Sys.Date(), units = "years")
mean(ages)
sd(ages)

# looking at individuals with thyroid cancer diagnosis and calculating demographics
ids <- as.numeric(colnames(mutants)[as.numeric(colnames(mutants)) %in% tc$person_id])
aad <- age_at_diag[age_at_diag$person_id %in% ids]
mean(aad$age_at_diagnosis)
sd(aad$age_at_diagnosis)
