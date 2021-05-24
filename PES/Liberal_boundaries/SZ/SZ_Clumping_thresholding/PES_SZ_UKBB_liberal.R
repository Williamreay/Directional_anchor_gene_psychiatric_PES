#############################################

## Association of PES and PRS with SZ in the UKBB cohort

## SZ cases and double the number of controls (randomly selected from MHQ participants)

## Liberal boundaries

## William Reay (2021)

#############################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/SZ_UKBB_liberal/")

library(data.table)
library(dplyr)
library(rcompanion)
library(corrplot)
library(ggplot2)
library(RColorBrewer)



## Read in concatenated PES file

All_SZ_scores <- fread("All_SZ_C_T_liberal_UKBB_PES.txt", header = T, sep = " ")

## Read in g.w. PRS

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores

PRS <- fread("../../PES/UKBB_SZ/SZ_UKBB_SZ_PGC3_PRS.all_score", header = T, sep = " ")

All_SZ_scores <- merge(All_SZ_scores, PRS, by = "IID")

## Read in phenotype and covariate data

SZ_pheno <- fread("../../UKBB_phenotypes/SZ_phenotypes/SZ_UKBB_pheno.txt", header = T)

SZ_cov <- fread("../../UKBB_phenotypes/SZ_phenotypes/Covariates_SZ_UKBB.txt", header = T)

## Merge all three df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                        list(All_SZ_scores, SZ_pheno, SZ_cov))

SZ_merged$Batch <- as.factor(SZ_merged$Batch)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

SZ_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 38:41)] <- lapply(SZ_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 38:41)], function(x) c(scale(x)))

## Define function to test each score by itself and extract beta, SE, and P

## List of sets to test

Scores_to_test <- as.list(colnames(SZ_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 38:41)]))

## Function

PES_PRS_SZ_UKBB <- function(v){
  Score_model <- glm(glue::glue('SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + {v} + Batch'), family = "binomial", data = SZ_merged)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

SZ_score <- sapply(Scores_to_test, PES_PRS_SZ_UKBB)

## Extract beta, se, z, and p value for each gene

SZ_extract <- apply(SZ_score, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

SZ_results <- data.frame()

for (i in 1:length(SZ_extract)) {
  SZ_results <- rbind(SZ_results, SZ_extract[[i]])
}
rownames(SZ_results) <- Scores_to_test

write.csv(SZ_results, file="Results/Liberal_SZ_UKBB_PES_PRS_results.csv", quote = F, row.names = T)


## Get variance explained on the liability scale (0.7% prevalence)

SZ_PES_PRS_r2 <- function(v, k, p){
  Null_model <- glm(SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch, family = "binomial", data = SZ_merged)
  Score_model <- glm(glue::glue('SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + {v}'), family = "binomial", data = SZ_merged)
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

SZ_r2 <- sapply(Scores_to_test, SZ_PES_PRS_r2, k=0.007, p=0.33)

R2_liability_results <- data.frame()

for (i in 1:length(SZ_r2)) {
  R2_liability_results <- rbind(R2_liability_results, SZ_r2[[i]])
}
rownames(R2_liability_results) <- Scores_to_test

write.csv(R2_liability_results, file="Results/Liberal_SZ_Liability_scale_r2_PES_PRS.csv", row.names = T, quote = F)

