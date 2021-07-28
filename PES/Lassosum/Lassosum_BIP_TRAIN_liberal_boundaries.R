#############################################

## Association of PES and PRS with BIP in the UKBB cohort - TRAIN

## Penalised regression

## Liberal boundaries

## BIP cases and double the number of controls (randomly selected from MHQ participants)

## William Reay (2021)

#############################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/BIP_TRAIN_lassosum_liberal/")

library(data.table)
library(dplyr)
library(rcompanion)
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(viridis)


## Read in concatenated PES/PRS file

All_BIP_scores <- fread("Liberal_all_BIP_TRAIN_UKBB_lassosum.txt", header = T, sep = " ")

Colnames_all_BIP_scores <- make.unique(names(All_BIP_scores))

colnames(All_BIP_scores) <- Colnames_all_BIP_scores

## Read in phenotype and covariate data

BIP_pheno <- fread("../../UKBB_phenotypes/BIP_phenotypes/BIP_TRAINING_UKBB_pheno.txt", header = T)

BIP_cov <- fread("../../UKBB_phenotypes/BIP_phenotypes/Covariates_TRAIN_BIP_UKBB.txt", header = T)

## Merge all three df

BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                     list(All_BIP_scores, BIP_pheno, BIP_cov))

BIP_merged$Batch <- as.factor(BIP_merged$Batch)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)] <- lapply(BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)], function(x) c(scale(x)))

## Define function to test each score by itself and extract beta, SE, and P

## List of sets to test

Scores_to_test <- as.list(colnames(BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)]))

## Function

PES_PRS_BIP_UKBB <- function(v){
  Score_model <- glm(glue::glue('BIP ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + {v} + Batch'), family = "binomial", data = BIP_merged)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

BIP_score <- sapply(Scores_to_test, PES_PRS_BIP_UKBB)

## Extract beta, se, z, and p value for each gene

BIP_extract <- apply(BIP_score, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

BIP_results <- data.frame()

for (i in 1:length(BIP_extract)) {
  BIP_results <- rbind(BIP_results, BIP_extract[[i]])
}
rownames(BIP_results) <- Scores_to_test

write.csv(BIP_results, file="Results/Liberal_Lassosum_BIP_UKBB_PES_PRS_results.csv", quote = F, row.names = T)


## Get variance explained on the liability scale (1% prevalence)

BIP_PES_PRS_r2 <- function(v, k, p){
  Null_model <- glm(BIP ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch, family = "binomial", data = BIP_merged)
  Score_model <- glm(glue::glue('BIP ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + {v}'), family = "binomial", data = BIP_merged)
  R2 <- nagelkerke(Score_model, null = Null_model)
  ## c.f. https://github.com/kn3in/ABC/blob/master/functions.R
  x <- qnorm(1 - k)
  z <- dnorm(x)
  i <- z / k
  cc <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
  theta <- i * ((p - k)/(1 - k)) * (i * ((p - k) / ( 1 - k)) - x)
  e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
  R2_liab <- cc * e * R2$Pseudo.R.squared.for.model.vs.null[3, ] / (1 + cc * e* theta * R2$Pseudo.R.squared.for.model.vs.null[3, ])
  return(R2_liab)
}

## Iterate function over list of genes

BIP_r2 <- sapply(Scores_to_test, BIP_PES_PRS_r2, k=0.01, p=0.33)

R2_liability_results <- data.frame()

for (i in 1:length(BIP_r2)) {
  R2_liability_results <- rbind(R2_liability_results, BIP_r2[[i]])
}
rownames(R2_liability_results) <- Scores_to_test

write.csv(R2_liability_results, file="Results/Liberal_Lassosum_BIP_Liability_scale_r2_PES_PRS.csv", row.names = T, quote = F)
