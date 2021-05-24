#############################################

## Association of PES and PRS with SZ in the ASRB cohort (sumstats with ASRB removed)

## Replicating scores from UKBB using PGC3 EUR weights with ASRB removed

## William Reay (2021)

#############################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/Replication_scores/ASRB_replication_full_cohort_ex_ASRB/")

library(data.table)
library(readxl)
library(dplyr)
library(rcompanion)

## Read in scores

All_SZ_scores <- fread("ASRB_replication_PES_PRS_combined.txt", header = T)

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores


## Read in phenotype and covariate data

SZ_pheno <- read_excel("ASRB_case_vs_control.xlsx")

## Merge case vs control df

ASRB_case_control <- merge(All_SZ_scores, SZ_pheno, by = "IID")


## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

ASRB_case_control[,c(6,12,18,24)] <- lapply(ASRB_case_control[,c(6,12,18,24)], function(x) c(scale(x)))


## Define function to test each score by itself and extract beta, SE, and P

## List of sets to test

Scores_to_test <- as.list(colnames(ASRB_case_control[,c(6,12,18,24)]))

## Function

PES_PRS_SZ_ASRB <- function(v){
  Score_model <- glm(glue::glue('PHENOTYPE ~ SEX + PC1 + PC2 + PC3 + {v}'), family = "binomial", data = ASRB_case_control)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

SZ_score <- sapply(Scores_to_test, PES_PRS_SZ_ASRB)

## Extract beta, se, z, and p value for each gene

SZ_extract <- apply(SZ_score, 2, function(x) return(as.data.frame(x$coefficients)[6,1:4]))

SZ_results <- data.frame()

for (i in 1:length(SZ_extract)) {
  SZ_results <- rbind(SZ_results, SZ_extract[[i]])
}
rownames(SZ_results) <- Scores_to_test

write.csv(SZ_results, file="SZ_ASRB_replication_PES_PRS_results.csv", quote = F, row.names = T)

## Liability scale

SZ_PES_PRS_r2 <- function(v, k, p){
  Null_model <- glm(PHENOTYPE ~ SEX + PC1 + PC2 + PC3, family = "binomial", data = ASRB_case_control)
  Score_model <- glm(glue::glue('PHENOTYPE ~ SEX + PC1 + PC2 + PC3 + {v}'), family = "binomial", data = ASRB_case_control)
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


SZ_r2 <- sapply(Scores_to_test, SZ_PES_PRS_r2, k=0.007, p=0.6286982)

R2_liability_results <- data.frame()

for (i in 1:length(SZ_r2)) {
  R2_liability_results <- rbind(R2_liability_results, SZ_r2[[i]])
}
rownames(R2_liability_results) <- Scores_to_test

write.csv(R2_liability_results, file="Results/ASRB_rep_Lassosum_ASRB_removed_SZ_Liability_scale_r2_PES_PRS.csv", row.names = T, quote = F)

## Test association with endophenotypes ##

## Merge with CD vs CS df

GoM <- read.csv("CD_vs_CS_GoM1.csv", header = T)

ASRB_GoM <- merge(All_SZ_scores, GoM, by = "IID")

ASRB_GoM[,c(6,12,18,24)] <- lapply(ASRB_GoM[,c(6,12,18,24)], function(x) c(scale(x)))

## Merge with clozapine df

Cloz <- read.table("Clozapine_ASRB_new_df.txt", header = T)

Cloz_merge <- merge(Cloz, SZ_pheno, by = "IID")

ASRB_cloz <- merge(All_SZ_scores, Cloz_merge, by = "IID")

ASRB_cloz[,c(6,12,18,24)] <- lapply(ASRB_cloz[,c(6,12,18,24)], function(x) c(scale(x)))

## GoM function

## Function

PES_PRS_SZ_GoM <- function(v){
  Score_model <- glm(glue::glue('GoM_coded ~ SEX + PC1 + PC2 + PC3 + {v}'), family = "binomial", data = ASRB_GoM)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

GoM_score <- sapply(Scores_to_test, PES_PRS_SZ_GoM)

## Extract beta, se, z, and p value for each gene

GoM_extract <- apply(GoM_score, 2, function(x) return(as.data.frame(x$coefficients)[6,1:4]))

GoM_results <- data.frame()

for (i in 1:length(GoM_extract)) {
  GoM_results <- rbind(GoM_results, GoM_extract[[i]])
}
rownames(GoM_results) <- Scores_to_test

write.csv(GoM_results, file="Results/SZ_ASRB_CD_CS_PES_PRS_results.csv", quote = F, row.names = T)

## Clozapine function

## Function

PES_PRS_SZ_Cloz <- function(v){
  Score_model <- glm(glue::glue('Clozapine ~ SEX + PC1 + PC2 + PC3 + {v}'), family = "binomial", data = ASRB_cloz)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

Cloz_score <- sapply(Scores_to_test, PES_PRS_SZ_Cloz)

## Extract beta, se, z, and p value for each gene

Cloz_extract <- apply(Cloz_score, 2, function(x) return(as.data.frame(x$coefficients)[6,1:4]))

Cloz_results <- data.frame()

for (i in 1:length(Cloz_extract)) {
  Cloz_results <- rbind(Cloz_results, Cloz_extract[[i]])
}
rownames(Cloz_results) <- Scores_to_test

write.csv(Cloz_results, file="Results/SZ_ASRB_clozapine_PES_PRS_results.csv", quote = F, row.names = T)
