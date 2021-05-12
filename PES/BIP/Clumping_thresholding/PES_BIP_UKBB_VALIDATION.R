#############################################

## Association of PES and PRS with BIP in the UKBB cohort - VALIDATION

## BIP cases and double the number of controls (randomly selected from MHQ participants)

## William Reay (2021)

#############################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/TEST_BIP/")

library(data.table)
library(dplyr)
library(rcompanion)
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)


## Read in concatenated PES/PRS file

All_BIP_scores <- fread("BIP_validation_PES_PRS_UKBB.txt", header = T, sep = " ")

Colnames_all_BIP_scores <- make.unique(names(All_BIP_scores))

colnames(All_BIP_scores) <- Colnames_all_BIP_scores

## Read in phenotype and covariate data

BIP_pheno <- fread("../../UKBB_phenotypes/BIP_phenotypes/BIP_TEST_UKBB_pheno.txt", header = T)

BIP_cov <- fread("../../UKBB_phenotypes/BIP_phenotypes/Covariates_TEST_BIP_UKBB.txt", header = T)

## Merge all three df

BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                     list(All_BIP_scores, BIP_pheno, BIP_cov))

BIP_merged$Batch <- as.factor(BIP_merged$Batch)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

BIP_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)] <- lapply(BIP_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)], function(x) c(scale(x)))

## Define function to test each score by itself and extract beta, SE, and P

## List of sets to test

Scores_to_test <- as.list(colnames(BIP_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)]))

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

write.csv(BIP_results, file="Results/VALIDATION_BIP_UKBB_PES_PRS_results.csv", quote = F, row.names = T)


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

write.csv(R2_liability_results, file="Results/VALIDATION_BIP_Liability_scale_r2_PES_PRS.csv", row.names = T, quote = F)

## Make forest plot of BIP OR per SD

BIP_results$Scores <- Scores_to_test


Plot_df <-  BIP_results %>% filter(Scores == "FES_0_005" |
                                     Scores == "RPS17_0_05" | Scores == "GRIN2A_0_005" |
                                     Scores == "CACNA1C_0_005" | Scores == "PCCB_0_05" |
                                     Scores == "FADS1_1")

Plot_df$Score_name <- c("CACNA1C netowrk PES", "FADS1 network PES", 
                        "FES network PES", "GRIN2A network PES",
                        "PCCB network PES", "RPS17 network PES")

## Plot results

Plot_df$OR <- exp(Plot_df$Estimate)
Plot_df$UOR <- Plot_df$OR + 1.96*Plot_df$`Std. Error`
Plot_df$LOR <- Plot_df$OR - 1.96*Plot_df$`Std. Error`

FP_BIP <- ggplot(data = Plot_df, aes(x=Score_name, y=OR, ymin=LOR, ymax=UOR, colour=Score_name)) +
  geom_pointrange() +
  geom_hline(yintercept=1, lty=2) +
  coord_flip() +
  ylab("Bipolar disorder OR [95% CI] per SD increase in PES") +
  theme_bw() +
  theme(legend.position = "null", axis.title.y = element_blank())

FP_BIP + scale_color_brewer(palette="Dark2")
