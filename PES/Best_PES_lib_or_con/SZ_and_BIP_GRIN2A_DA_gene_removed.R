############################

## GRIN2A lassosum PES - GRIN2A region removed

############################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/")

library(data.table)
library(dplyr)
library(rcompanion)
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)

## SZ ##

## Read in concatenated PES/PRS file

All_SZ_scores <- fread("../../Best_DA_gene_PES/DA_gene_removed/SZ_UKBB_LASSO_SZ_GRIN2A_DA_gene_removed.validate.results.txt", header = T, sep = " ")

## Read in phenotype and covariate data

SZ_pheno <- fread("../../UKBB_phenotypes/SZ_phenotypes/SZ_UKBB_pheno.txt", header = T)

SZ_cov <- fread("../../UKBB_phenotypes/SZ_phenotypes/Covariates_SZ_UKBB.txt", header = T)

## Merge all three df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_SZ_scores, SZ_pheno, SZ_cov))

SZ_merged$Batch <- as.factor(SZ_merged$Batch)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

SZ_merged[,5] <- lapply(SZ_merged[,5], function(x) c(scale(x)))

## Test associatiom

Score_model_SZ <- glm(SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + SZ_GRIN2A_DA_gene_removed + Batch, family = "binomial",
                      data = SZ_merged)

## BIP ##

All_BIP_scores <- fread("../../Best_DA_gene_PES/DA_gene_removed/BIP_TRAIN_UKBB_LASSO_BIP_GRIN2A_DA_gene_removed.validate.results.txt",
                        sep = " ", header = T)

BIP_pheno <- fread("../../UKBB_phenotypes/BIP_phenotypes/BIP_TRAINING_UKBB_pheno.txt", header = T)

BIP_cov <- fread("../../UKBB_phenotypes/BIP_phenotypes/Covariates_TRAIN_BIP_UKBB.txt", header = T)

BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_BIP_scores, BIP_pheno, BIP_cov))

BIP_merged$Batch <- as.factor(BIP_merged$Batch)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

BIP_merged[,5] <- lapply(BIP_merged[,5], function(x) c(scale(x)))

Score_model_BIP <- glm(BIP ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + BIP_DA_Gene_GRIN2A_removed + Batch, family = "binomial",
                      data = BIP_merged)
