########################################

## Testing the association between SZ and BIP PES with psychological indicies of a general population cohort

## Hunter community study

## William Reay (2021)

#####################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/")

## Load dependencies

library(dplyr)
library(data.table)
library(easyGgplot2)
library(ggpubr)
library(readxl)

## Load HCS data and identify individuals with non-missing K10 psychological distress, CESD scale, and self-reported depression/anxiety

HCS_pheno <- read.csv("../Hunter_cohort/Phenotype_data/mcv1.csv", header = T)

HCS_pheno <- HCS_pheno %>% filter(!is.na(k10w1) & !is.na(CESDw1) & !is.na(Depression_Anxietyw1))

## Load covariates

HCS_cov <- fread("../Hunter_cohort/Eigenvectors/HCS_PC1_PC5.tab.txt", header = F)

HCS_cov <- rename(HCS_cov, "ID"="V1", "PC1"="V3", "PC2"="V4", "PC3"="V5",
                  "PC4"="V6", "PC5"="V7")

Electoral_HCSID_conversion <- read_excel("~/cloudstor//Pneumonia_cytokine_lung_function/HCS_PES_profiles/Electoral_HCSID_conversion.xlsx")

## Derive EIDs from Hxxxx IDs

HCS_cov <- merge(HCS_cov, Electoral_HCSID_conversion, by = "ID")

HCS_cov <- rename(HCS_cov, "IID"="electoralId")


## Load SZ scores and append SZ prefix to them

SZ_scores <- fread("Best_DA_gene_PES/HCS_cohort_profiling/SZ_HCS_scores/Combined_SZ_HCS_scores.txt", header = T)

Colnames_all_SZ_scores <- make.unique(names(SZ_scores))

colnames(SZ_scores) <- Colnames_all_SZ_scores

SZ_scores <- SZ_scores %>% rename_all( ~ paste0("SZ_", .x))
SZ_scores <- rename(SZ_scores, "IID"="SZ_IID")

## Load BIP scores and append BIP prefix to them

BIP_scores <- fread("Best_DA_gene_PES/HCS_cohort_profiling/BIP_HCS_scores/Combined_BIP_HCS_scores.txt", header = T)

Colnames_all_BIP_scores <- make.unique(names(BIP_scores))

colnames(BIP_scores) <- Colnames_all_BIP_scores

BIP_scores <- BIP_scores %>% rename_all( ~ paste0("BIP_", .x))
BIP_scores <- rename(BIP_scores, "IID"="BIP_IID")

## Merge all three

HCS_pheno <- rename(HCS_pheno, "IID"="electoralId")

SZ_and_BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                            list(SZ_scores, BIP_scores, HCS_pheno, HCS_cov))

SZ_and_BIP_merged <- rename(SZ_and_BIP_merged, "SZ_RPS17_PES_AVG"="SZ_RPS17$best.beta_AVG")

## Scale the BIP and SZ PES/PRS

SZ_and_BIP_merged[,c(5, 15, 25, 35, 45, 55, 65,74, 79, 84, 89, 94, 99, 104)] <- lapply(SZ_and_BIP_merged[,c(5, 15, 25, 35, 45, 55, 65,74, 79, 84, 89, 94, 99, 104)], function(x) c(scale(x)))

## Define function to test each score by itself and extract beta, SE, and P

## List of sets to test

Scores_to_test <- as.list(colnames(SZ_and_BIP_merged[,c(5, 15, 25, 35, 45, 55, 65, 74, 79, 84, 89, 94, 99, 104)]))

## Make psych distress cut-off of 17

SZ_and_BIP_merged$K10_cutoff <- as.factor(ifelse(SZ_and_BIP_merged$k10w1 >= 17, 1, 0))

## Function

HCS_k10 <- function(v){
  Score_model <- glm(glue::glue('K10_cutoff ~ sex.x + bage + PC1 + PC2 + PC3 + PC4 + PC5 + {v}'), family = "binomial", data = SZ_and_BIP_merged)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

K10_score <- sapply(Scores_to_test, HCS_k10)

## Extract beta, se, z, and p value for each gene

K10_extract <- apply(K10_score, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

K10_results <- data.frame()

for (i in 1:length(K10_extract)) {
  K10_results <- rbind(K10_results, K10_extract[[i]])
}
rownames(K10_results) <- Scores_to_test

write.csv(K10_results, file="Best_DA_gene_PES/HCS_cohort_profiling/Association_results/K10_PES_PRS_BIP_and_SZ.csv", quote = F, row.names = T)

## Look at self-reported anxiety/depression

HCS_anx_dep <- function(v){
  Score_model <- glm(glue::glue('Depression_Anxietyw1 ~ sex.x + bage + PC1 + PC2 + PC3 + PC4 + PC5 + {v}'), family = "binomial", data = SZ_and_BIP_merged)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

Anx_score <- sapply(Scores_to_test, HCS_anx_dep)

## Extract beta, se, z, and p value for each gene

Anx_extract <- apply(Anx_score, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

Anx_results <- data.frame()

for (i in 1:length(Anx_extract)) {
  Anx_results <- rbind(Anx_results, Anx_extract[[i]])
}
rownames(Anx_results) <- Scores_to_test

write.csv(Anx_results, file="Best_DA_gene_PES/HCS_cohort_profiling/Association_results/Self_report_anxiety_depression_PES_PRS_BIP_and_SZ.csv", quote = F, row.names = T)

## Test FADS1 SZ PES assoc for PRS covariation

FADS1_anx_PRS <- glm(Depression_Anxietyw1 ~ sex.x + bage + PC1 + PC2 + PC3 + PC4 + PC5 + SZ_SZ_PRS_ALL_AVG +
                       SZ_FADS1_PES_ALL_AVG, family = "binomial", data = SZ_and_BIP_merged)

## Look at CES-D > 20

SZ_and_BIP_merged$CESD_cutoff <- as.factor(ifelse(SZ_and_BIP_merged$CESDw1 >= 20, 1, 0))

HCS_cesd <- function(v){
  Score_model <- glm(glue::glue('CESD_cutoff ~ sex.x + bage + PC1 + PC2 + PC3 + PC4 + PC5 + {v}'), family = "binomial", data = SZ_and_BIP_merged)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

cesd_score <- sapply(Scores_to_test, HCS_cesd)

## Extract beta, se, z, and p value for each gene

cesd_extract <- apply(cesd_score, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

cesd_results <- data.frame()

for (i in 1:length(cesd_extract)) {
  cesd_results <- rbind(cesd_results, cesd_extract[[i]])
}
rownames(cesd_results) <- Scores_to_test

write.csv(cesd_results, file="Best_DA_gene_PES/HCS_cohort_profiling/Association_results/CESD_PES_PRS_BIP_and_SZ.csv", quote = F, row.names = T)

