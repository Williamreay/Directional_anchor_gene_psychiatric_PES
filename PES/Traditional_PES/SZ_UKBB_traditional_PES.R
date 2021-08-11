#############################################

## Traditional PES in the UKBB cohort - SZ

## SZ cases and double the number of controls (randomly selected from MHQ participants)

## William Reay (2021)

#############################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Updated_SZ_BIP_PES_pathways/SZ_TRAINING_TRADITIONAL_PES/")

library(data.table)
library(dplyr)
library(rcompanion)
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)


## Read in concatenated PES file

All_SZ_scores <- fread("All_SZ_traditional_PES_UKBB.txt", header = T, sep = " ")

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores

## Read in phenotype and covariate data

SZ_pheno <- fread("../../UKBB_phenotypes/SZ_phenotypes/SZ_UKBB_pheno.txt", header = T)

SZ_cov <- fread("../../UKBB_phenotypes/SZ_phenotypes/Covariates_SZ_UKBB.txt", header = T)

## Read in PRS

PRS <- fread("../../Best_DA_gene_PES/SZ_best_con_or_lib/SZ_UKBB_TRAIN_score_files/SZ_UKBB_LASSO_SZ_PGC3_PRS.validate.results.txt", header = T)

## Merge all four df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_SZ_scores, PRS, SZ_pheno, SZ_cov))

SZ_merged$Batch <- as.factor(SZ_merged$Batch)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

SZ_merged[,c(3,7,10,14,18,21,25,29,32,36,41)] <- lapply(SZ_merged[,c(3,7,10,14,18,21,25,29,32,36,41)], function(x) c(scale(x)))

## List of sets to test

Scores_to_test <- as.list(colnames(SZ_merged[,c(3,7,10,14,18,21,25,29,32,36,41)]))

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

write.csv(SZ_results, file="Results/SZ_UKBB_PES_PRS_results.csv", quote = F, row.names = T)

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

write.csv(R2_liability_results, file="Results/SZ_Liability_scale_r2_PES_PRS.csv", row.names = T, quote = F)

## Test covariation of PRS

## Spermatogenesis

Spermatogenesis <- glm(SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS + SPERMATOGENESIS_0_005,
                       family="binomial", data = SZ_merged)

anova(Spermatogenesis, test="Chisq")

## Neuronal system

Neuronal_system <- glm(SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS + Neuronal_system_1,
                       family="binomial", data = SZ_merged)

anova(Neuronal_system, test="Chisq")

## Insulin secretion

Insulin <- glm(SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS + Insulin_secretion_0_5,
                       family="binomial", data = SZ_merged)

anova(Insulin, test="Chisq")

## Correlation plot

Cor_SZ <- cor(SZ_merged[,c(3,7,10,14,18,21,25,29,32,36,41)])

pval <- psych::corr.test(SZ_merged[,c(3,7,10,14,18,21,25,29,32,36,41)], adjust = "none")$p

Sig_corr_plot <- corrplot(Cor_SZ, tl.cex = 0.5, p.mat = pval, insig = "blank",
                          sig.level = 0.05, type="upper", method="number",
                          tl.col="black", tl.srt=45)

SZ_results$Scores <- Scores_to_test

SZ_results$Score_name <- c("Nuclear beta catenin signalling", "Cardiac muscle secretion", 
                        "CK1 pathway", "CRMPs in SEMA3A signalling",
                        "Hypertrophic cardiomyopathy", "Regulation of insulin secretion",
                        "mRNA stability", "Neuronal system", "PIP3 signalling in cardiac myocytes",
                        "Spermatogenesis", "PRS")

Plot_df <- SZ_results

Plot_df <- Plot_df %>% filter(Score_name != "PRS")

Plot_df$OR <- exp(Plot_df$Estimate)
Plot_df$UOR <- Plot_df$OR + 1.96*Plot_df$`Std. Error`
Plot_df$LOR <- Plot_df$OR - 1.96*Plot_df$`Std. Error`

FP_SZ <- ggplot(data = Plot_df, aes(x=Score_name, y=OR, ymin=LOR, ymax=UOR, colour=Score_name)) +
  geom_pointrange() +
  geom_hline(yintercept=1, lty=2) +
  coord_flip() +
  ylab("Schizophrenia OR [95% CI] per SD increase in PES") +
  theme_bw() +
  theme(legend.position = "null", axis.title.y = element_blank())
