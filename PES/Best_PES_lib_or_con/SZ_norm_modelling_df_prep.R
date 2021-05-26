###############################

## Preparing SZ df for normative modelling

##############################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/UKBB_SZ/")

library(dplyr)
library(data.table)
library(quantreg)
library(easyGgplot2)

set.seed(57373)

All_SZ_scores <- fread("Combined_SZ_all_PES_PRS_UKBB.txt", header = T, sep = " ")

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores

## Read in phenotype and covariate data

SZ_pheno <- fread("../../UKBB_phenotypes/SZ_phenotypes/SZ_UKBB_pheno.txt", header = T)

SZ_cov <- fread("../../UKBB_phenotypes/SZ_phenotypes/Covariates_SZ_UKBB.txt", header = T)

## Merge all three df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_SZ_scores, SZ_pheno, SZ_cov))

SZ_merged$Phenotype <- ifelse(SZ_merged$SZ == 1, "SZ", "HC")

SZ_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)] <- lapply(SZ_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)], function(x) c(scale(x)))

SZ_merged$Batch <- as.factor(SZ_merged$Batch)


## Quantile regression

Controls <- SZ_merged %>% filter(SZ == 0)

Cases <- SZ_merged %>% filter(SZ == 1)

Multi_FES <- rq(FES_0_5 ~ PRS_0_5, data = Controls, tau = c(0.1, 0.5, 0.9))

Multi_sum <- summary(Multi_FES, se="boot")

colors <- c("firebrick1", "royalblue", "goldenrod")
plot(FES_0_5 ~ PRS_0_5, data = SZ_merged, pch = 16, cex = 0.7, main = "FES network PES ~ PRS")
for (j in 1:ncol(Multi_FES$coefficients)) {
  abline(coef(Multi_FES)[, j], col = colors[j])
}

## Predict in cases

SZ_pred <- predict(Multi_FES, data.frame(Cases), interval = "prediction", level = 0.1)

Cases$Tau_0_1_FES <- as.data.frame(SZ_pred[, 1])

Cases$Tau_0_5_FES <- as.data.frame(SZ_pred[, 2])

Cases$Tau_0_9_FES <- as.data.frame(SZ_pred[, 3])

## Get SD in controls

sd(Controls$FES_0_5)

Cases$Norm_0_9 <- (Cases$FES_0_5 - Cases$Tau_0_9_FES)/0.000323711

Cases$Norm_0_5 <- (Cases$FES_0_5 - Cases$Tau_0_5_FES)/0.000323711

Cases$Norm_0_1 <- (Cases$FES_0_5 - Cases$Tau_0_1_FES)/0.000323711

