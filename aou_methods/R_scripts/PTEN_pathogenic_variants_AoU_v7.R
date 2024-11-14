install.packages('R.utils')
install.packages("eeptools")
library(data.table)
library(eeptools)


all <- fread("./PHENOTYPES/srWGS_all_demographics.csv")
tc <- fread("./PHENOTYPES/thyroid_cancer_srWGS.csv") 

ret = as.matrix(fread("./GENOTYPES/PTEN_genotypes.tsv.gz"))
dim(ret)

# last column has NAs removing

ret <- ret[, 1:(ncol(ret)-1)]

# getting the list of pathogenic variants from ClinVar table
vars <- fread("./GENES/PTEN.tsv.gz")

# pathalleles <- c(paste0("chr", vars$Chromosome[1], ":", vars$Start, ":", vars$ReferenceAlleleVCF, ":", vars$AlternateAlleleVCF),
#                 paste0("chr", vars$Chromosome[1], ":", vars$Start, ":", vars$AlternateAlleleVCF, ":", vars$ReferenceAlleleVCF))

# manually curated variants against the pull for DICER1 gene from clinVar

# ret[, "variant"][!(ret[, "variant"] %in% pathalleles)]
# keep <- ret[, "variant"][ret[, "variant"] %in% pathalleles]




# manually curated list of P/LP variants
keep <- c("chr10:87864334:G:C", "chr10:87864506:A:T", "chr10:87864514:A:T", "chr10:87864546:C:T", "chr10:87894049:T:G", "chr10:87925550:T:C",
  "chr10:87931090:G:A", "chr10:87933061:T:C", "chr10:87933073:G:A", "chr10:87933147:C:T", "chr10:87933148:G:A", "chr10:87933165:T:C",
  "chr10:87933222:T:C", "chr10:87952142:C:T", "chr10:87952207:G:GT", "chr10:87957915:C:T", "chr10:87960892:A:G", "chr10:87961095:C:T",
  "chr10:87961100:CT:C", "chr10:87965292:G:GAC")
 
# filtering for P/LP variants that matched
ret <- ret[ret[, "variant"] %in% keep, ]

# removing non-genotype column
varnames <- ret[, "variant"]

ret <- ret[,!colnames(ret) %in% c("variant")]

rownames(ret) <- varnames

# identifying individuals with P/LP variants
muts <- data.frame(variants = varnames,
                   counts = numeric(length(varnames)))

sum(is.na(colnames(ret)))

mutants <- ret[, ret[1,] != "0/0" & ret[1,] != "./.", drop = F]
muts$counts[1] <- dim(mutants)[2]

for (i in 2:nrow(ret)){

  tmp <- ret[, ret[i,] != "0/0" & ret[i,] != "./.", drop = F]
  mutants <- cbind(mutants, tmp)
  muts$counts[i] <- dim(tmp)[2]
}

print(muts)


write.csv(mutants, "./GENOTYPES/patients_with_pathogenic_mutations_in_PTEN.csv")
write.csv(muts, "./GENOTYPES/COUNTS_of_pathogenic_mutations_in_PTEN.csv")


as.numeric(colnames(mutants)) %in% clinical$person_id
as.numeric(colnames(mutants)) %in% tc$person_id


# demographics for all APC mutants
all_PTEN <- all[all$person_id %in% as.numeric(colnames(mutants)), ]
ages <- age_calc(as.Date(all_PTEN$date_of_birth), enddate = Sys.Date(), units = "years")
mean(ages)
sd(ages)
table(all_PTEN$sex_at_birth)

# looking at individuals with thyroid cancer diagnosis and calculating demographics
ids <- as.numeric(colnames(mutants)[as.numeric(colnames(mutants)) %in% tc$person_id])
aad <- age_at_diag[age_at_diag$person_id %in% ids]
mean(aad$age_at_diagnosis)
sd(aad$age_at_diagnosis)

