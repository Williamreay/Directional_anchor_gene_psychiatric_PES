#############################################

## GRIN2A DA-gene removed

## SZ/BIP cases and double the number of controls (randomly selected from MHQ participants)

## William Reay (2021)

#############################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/")

library(data.table)
library(dplyr)
library(rcompanion)
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)

## SZ

## Read in phenotype and covariate data

SZ_pheno <- fread("../../UKBB_phenotypes/SZ_phenotypes/SZ_UKBB_pheno.txt", header = T)

SZ_cov <- fread("../../UKBB_phenotypes/SZ_phenotypes/Covariates_SZ_UKBB.txt", header = T)

SZ_GRIN2A <- fread("../../Best_DA_gene_PES/DA_gene_removed/SZ_UKBB_LASSO_SZ_GRIN2A_DA_gene_removed.validate.results.txt")

## Merge all three df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(SZ_GRIN2A, SZ_pheno, SZ_cov))

SZ_merged$Batch <- as.factor(SZ_merged$Batch)

SZ_merged$SZ_GRIN2A_DA_gene_removed <- as.numeric(scale(SZ_merged$SZ_GRIN2A_DA_gene_removed))

## Model

GRIN2A_SZ_mod <- glm(SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + SZ_GRIN2A_DA_gene_removed + Batch,
                     family = "binomial", data = SZ_merged)

## BIP

BIP_GRIN2A <- fread("../../Best_DA_gene_PES/DA_gene_removed/BIP_TRAIN_UKBB_LASSO_BIP_GRIN2A_DA_gene_removed.validate.results.txt", header = T)

BIP_pheno <- fread("../../UKBB_phenotypes/BIP_phenotypes/BIP_TRAINING_UKBB_pheno.txt", header = T)

BIP_cov <- fread("../../UKBB_phenotypes/BIP_phenotypes/Covariates_TRAIN_BIP_UKBB.txt", header = T)

## Merge all three df

BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                     list(BIP_GRIN2A, BIP_pheno, BIP_cov))

BIP_merged$Batch <- as.factor(BIP_merged$Batch)

BIP_merged$GRIN2A_BIP_removed_DA <- as.numeric(scale(BIP_merged$GRIN2A_BIP_removed_DA))

GRIN2A_BIP_mod <- glm(BIP ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + GRIN2A_BIP_removed_DA + Batch,
                     family = "binomial", data = BIP_merged)
