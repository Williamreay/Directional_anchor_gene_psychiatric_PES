######################################

## Replicating PES/PRS in the HCS cohort

## Lipid correlations 

## William Reay (2021)

######################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem/")

library(dplyr)
library(data.table)
library(readxl)

## Read in PES/PRS

## Load SZ scores and append SZ prefix to them

SZ_scores <- fread("../HCS_cohort_profiling/SZ_HCS_scores/Combined_SZ_HCS_scores.txt", header = T)

Colnames_all_SZ_scores <- make.unique(names(SZ_scores))

colnames(SZ_scores) <- Colnames_all_SZ_scores

SZ_scores <- SZ_scores %>% rename_all( ~ paste0("SZ_", .x))
SZ_scores <- rename(SZ_scores, "IID"="SZ_IID")

## Load BIP scores and append BIP prefix to them

BIP_scores <- fread("../HCS_cohort_profiling/BIP_HCS_scores/Combined_BIP_HCS_scores.txt", header = T)

Colnames_all_BIP_scores <- make.unique(names(BIP_scores))

colnames(BIP_scores) <- Colnames_all_BIP_scores

BIP_scores <- BIP_scores %>% rename_all( ~ paste0("BIP_", .x))
BIP_scores <- rename(BIP_scores, "IID"="BIP_IID")


## Read in HCS data and retain those indidviduals with non MH conditions

HCS_pheno <- read.csv("~/Desktop/Hunter_cohort/Phenotype_data/mcv1.csv", header = T)

HCS_pheno <- HCS_pheno %>% filter(!is.na(k10w1) & !is.na(CESDw1) & !is.na(Depression_Anxietyw1))

## Load covariates

HCS_cov <- fread("~/Desktop/Hunter_cohort/Eigenvectors/HCS_PC1_PC5.tab.txt", header = F)

HCS_cov <- rename(HCS_cov, "ID"="V1", "PC1"="V3", "PC2"="V4", "PC3"="V5",
                  "PC4"="V6", "PC5"="V7")

Electoral_HCSID_conversion <- read_excel("~/cloudstor/Pneumonia_cytokine_lung_function/HCS_PES_profiles/Electoral_HCSID_conversion.xlsx")

HCS_cov <- merge(HCS_cov, Electoral_HCSID_conversion, by = "ID")

HCS_cov <- rename(HCS_cov, "IID"="electoralId")

HCS_pheno <- rename(HCS_pheno, "IID"="electoralId")

SZ_and_BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                            list(SZ_scores, BIP_scores, HCS_pheno, HCS_cov))

## Test association with lipids, retain individuals with non-missing HDL, LDL, and TG

SZ_and_BIP_merged <- SZ_and_BIP_merged %>% filter(!is.na(HDLT) & !is.na(LDLT) & !is.na(TRIG))

## Scale FADS1 PES and PRS

SZ_and_BIP_merged$SZ_FADS1_PES_ALL_AVG <- as.numeric(scale(SZ_and_BIP_merged$SZ_FADS1_PES_ALL_AVG))
SZ_and_BIP_merged$BIP_FADS1_PES_AVG <- as.numeric(scale(SZ_and_BIP_merged$BIP_FADS1_PES_AVG))
SZ_and_BIP_merged$SZ_SZ_PRS_ALL_AVG <- as.numeric(scale(SZ_and_BIP_merged$SZ_SZ_PRS_ALL_AVG))
SZ_and_BIP_merged$BIP_BIP_PRS_AVG <- as.numeric(scale(SZ_and_BIP_merged$BIP_BIP_PRS_AVG))


List_of_scores <- as.list(c("SZ_FADS1_PES_ALL_AVG", "BIP_FADS1_PES_AVG", "SZ_SZ_PRS_ALL_AVG", "BIP_BIP_PRS_AVG"))

## Test association with TG (mmol/L)

HCS_TG <- function(v){
  Score_model <- lm(glue::glue('TRIG ~ sex.x + bage + PC1 + PC2 + PC3 + PC4 + PC5 + {v}'), data = SZ_and_BIP_merged)
  return(summary(Score_model)) 
}

TG_test <- sapply(List_of_scores, HCS_TG)

TG_extract <- apply(TG_test, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

## Test association with HDL (mmol/L)

HCS_HDL <- function(v){
  Score_model <- lm(glue::glue('HDLT ~ sex.x + bage + PC1 + PC2 + PC3 + PC4 + PC5 + {v}'), data = SZ_and_BIP_merged)
  return(summary(Score_model)) 
}

HDL_test <- sapply(List_of_scores, HCS_HDL)

HDL_extract <- apply(HDL_test, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

## Test association with LDL (mmol/L)

HCS_LDL <- function(v){
  Score_model <- lm(glue::glue('LDLT ~ sex.x + bage + PC1 + PC2 + PC3 + PC4 + PC5 + {v}'), data = SZ_and_BIP_merged)
  return(summary(Score_model)) 
}

LDL_test <- sapply(List_of_scores, HCS_LDL)

LDL_extract <- apply(LDL_test, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))
