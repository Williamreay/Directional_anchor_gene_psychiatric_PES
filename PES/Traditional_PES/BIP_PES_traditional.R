#############################################

## Association of PES and PRS with BIP in the UKBB cohort - TRAIN

## BIP cases and double the number of controls (randomly selected from MHQ participants)

## William Reay (2021)

#############################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Updated_SZ_BIP_PES_pathways/BIP_TRAINING_TRADITIONAL_PES/")

library(data.table)
library(dplyr)
library(rcompanion)
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)

## Read in concatenated PES/PRS file

All_BIP_scores <- fread("All_BIP_trad_PES.txt", header = T, sep = " ")

Colnames_all_BIP_scores <- make.unique(names(All_BIP_scores))

colnames(All_BIP_scores) <- Colnames_all_BIP_scores

All_BIP_scores <- All_BIP_scores %>% select(IID, ECM_BIP_0_005_PES,
                                            EGFR_1_PES, FC_EPI_BIP_0_5_PES,
                                            IL8_CXCR2_0_005_PES, INSULIN_ACETYLCHOLINE_0_005_PES,
                                            INOSITOL_PHOSPHATE_0_005_PES, Ligand_gated_ion_0_5_PES,
                                            Lung_cancer_0_5_PES, MAPK_TRK_1_BIP_PES, Neuronal_system_BIP_PES,
                                            BIP_OPIOID_0_005_PES, Phosphatdylinotosol_0_005_PES,
                                            Progesterone_0_5_PES, Steroid_0_005_PES,
                                            TCR_RAS_0_5_PES)

## PRS

PRS <- fread("../../Best_DA_gene_PES/BIP_best_con_or_lib/BIP_UKBB_TRAIN_score_files/BIP_PRS_penalised_regression_train_score.txt", header = T)

## Read in phenotype and covariate data

BIP_pheno <- fread("../../UKBB_phenotypes/BIP_phenotypes/BIP_TRAINING_UKBB_pheno.txt", header = T)

BIP_cov <- fread("../../UKBB_phenotypes/BIP_phenotypes/Covariates_TRAIN_BIP_UKBB.txt", header = T)

## Merge all three df

BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                     list(All_BIP_scores, PRS, BIP_pheno, BIP_cov))

BIP_merged$Batch <- as.factor(BIP_merged$Batch)

BIP_merged[,c(2:17)] <- lapply(BIP_merged[,c(2:17)], function(x) c(scale(x)))

Scores_to_test <- as.list(colnames(BIP_merged[,c(2:17)]))

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

write.csv(BIP_results, file="Results/Traditional_PES_PRS_BIP.csv", quote = F, row.names = T)


## Get variance explained on the liability scale (0.7% prevalence)

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

write.csv(R2_liability_results, file="Results/Traditional_PES_PRS_BIP_Liability_scale_r2_PES_PRS.csv", row.names = T, quote = F)

## Ligand gated PES covariation

Ligand <- glm(BIP ~ Sex + Age + PC1 + PC2 + PC3 +
                PC4 + PC5 + Batch + PRS + Ligand_gated_ion_0_5_PES, family = "binomial", data = BIP_merged)

anova(Ligand, test = "Chisq")


## Neuronal PES covariation

Neuronal <- glm(BIP ~ Sex + Age + PC1 + PC2 + PC3 +
                PC4 + PC5 + Batch + PRS + Neuronal_system_BIP_PES, family = "binomial", data = BIP_merged)

anova(Neuronal, test = "Chisq")

