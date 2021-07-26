#############################

## Clustering the biochemical measures within the SZ and BIP cohorts

## William Reay (2021)

############################

library(dplyr)
library(data.table)
library(mclust)
library(tidyverse)

## Read in biochem df

Raw_biochem <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/UKBB_biochemical_pheWAS_IDs/Meausred_biochemical_data_UKBB.txt", header = T)

## Read in list of columns to retain (oestradiol excluded)

List_col <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem/Clustering_list.txt", header = F)

List_col <- as.list(List_col)

DF_biochem_raw <- as.data.frame(Raw_biochem)

df.subset <- DF_biochem_raw[, List_col[[1]]]

## Merge with the SZ cohort

SZ_pheno <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/SZ_phenotypes/All_SZ_unrelated_EUR_genotype_avail_UKBB.txt",header = T)

SZ_raw_biochem_merge <- merge(SZ_pheno, df.subset, by  ="f.eid")

colSums(is.na(SZ_raw_biochem_merge))

## Make SZ cohort with no missing values, remove rheumatoid factor and urine microalbumin due to high missingness

No_missing_SZ <- SZ_raw_biochem_merge %>% select(!(f.30820.0.0)) %>% select(!(f.30500.0.0))

No_missing_SZ <- na.omit(No_missing_SZ)

No_missing_SZ[ ,3:33] <- lapply(No_missing_SZ[,3:33], function(x) c(scale(x)))

BIC_all <- mclustBIC(No_missing_SZ[, 3:33])

plot(BIC_all)

SZ_biochem_complete_rows <- Mclust(No_missing_SZ[, 3:33], x = BIC_all)

summary(SZ_biochem_complete_rows , parameters = T)

plot(SZ_biochem_complete_rows, what ="classification", dimens = c(1,2,3))
