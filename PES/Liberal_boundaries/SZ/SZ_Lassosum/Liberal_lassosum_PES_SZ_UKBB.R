#############################################

## Association of penalised regression PES and PRS with SZ in the UKBB cohort

## SZ cases and double the number of controls (randomly selected from MHQ participants)

## Liberal boundaries

## William Reay (2021)

#############################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/SZ_UKBB_lassosum_liberal/")

library(data.table)
library(dplyr)
library(rcompanion)
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)


## Read in concatenated PES file

All_SZ_scores <- fread("All_SZ_lassosum_liberal_PES.txt", header = T, sep = " ")

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores

## Read in PRS

PRS <- fread("../../PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_PRS.validate.results.txt", header = T, sep = " ")

All_SZ_scores <- merge(All_SZ_scores, PRS, by = "IID")

## Read in phenotype and covariate data

SZ_pheno <- fread("../../UKBB_phenotypes/SZ_phenotypes/SZ_UKBB_pheno.txt", header = T)

SZ_cov <- fread("../../UKBB_phenotypes/SZ_phenotypes/Covariates_SZ_UKBB.txt", header = T)

## Merge all three df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_SZ_scores, SZ_pheno, SZ_cov))

SZ_merged$Batch <- as.factor(SZ_merged$Batch)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

SZ_merged[,c(5, 10, 15, 20, 25, 30, 34)] <- lapply(SZ_merged[,c(5, 10, 15, 20, 25, 30, 34)], function(x) c(scale(x)))

## Define function to test each score by itself and extract beta, SE, and P

## List of sets to test

Scores_to_test <- as.list(colnames(SZ_merged[,c(5, 10, 15, 20, 25, 30, 34)]))

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

write.csv(SZ_results, file="Results/Liberal_Lassosum_SZ_UKBB_PES_PRS_results.csv", quote = F, row.names = T)


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

write.csv(R2_liability_results, file="Results/Liberal_Lassosum_SZ_Liability_scale_r2_PES_PRS.csv", row.names = T, quote = F)


## Test covaration of FES PES with PRS

FES_full <-  glm(SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS + FES_PES, family = "binomial", data = SZ_merged)

Chi_sq_FES <- anova(FES_full, test = "Chisq")

FES_corr <- lm(PRS ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + FES_PES, data = SZ_merged)


## Test covariation with GRIN2A PES and also correlation

GRIN2A_full <-  glm(SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS + GRIN2A_PES, family = "binomial", data = SZ_merged)

Chi_sq_GRIN2A <- anova(GRIN2A_full, test = "Chisq")

GRIN2A_corr <- lm(GRIN2A_PES ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS, data = SZ_merged)

## RPS17 correlation and covariation

RPS17_full <-  glm(SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS + RPS17_PES, family = "binomial", data = SZ_merged)

Chi_sq_RPS17 <- anova(RPS17_full, test = "Chisq")

RPS17_corr <- lm(RPS17_PES ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS, data = SZ_merged)


## Correlation plot

Cor_SZ <- cor(SZ_merged[,c(5, 10, 15, 20, 25, 30, 34)])

pval <- psych::corr.test(SZ_merged[,c(5, 10, 15, 20, 25, 30, 34)], adjust = "none")$p

Sig_corr_plot <- corrplot(Cor_SZ, tl.cex = 0.5, p.mat = pval, insig = "blank",
                          sig.level = 0.05, type="upper", method="number",
                          tl.col="black", tl.srt=45)

Bon_corr_plot <- corrplot(Cor_SZ, tl.cex = 0.5, p.mat = pval, insig = "blank",
                          sig.level = 0.001020408, type="upper", method="number",
                          tl.col="black", tl.srt=45)


## Make forest plot of effect size per SD

SZ_results$Scores <- unlist(Scores_to_test)

Plot_df <- SZ_results %>% filter(Scores == "CACNA1C_PES" | Scores == "FADS1_PES" | Scores == "GRIN2A_PES" |
                                 Scores == "FES_PES" | Scores == "PCCB_PES" | Scores == "RPS17_PES")

Plot_df$Score_name <- c("CACNA1C network PES", "FADS1 network PES", 
                        "FES network PES", "GRIN2A network PES",
                        "PCCB network PES", "RPS17 network PES")

## Plot results

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

FP_SZ + scale_color_brewer(palette="YlOrRd")

## Remove cases who also meet a BIP criteria

Overlap <- fread("../../UKBB_phenotypes/IDs_meeting_criteria_for_SZ_and_BIP.txt", header = T)

Overlap <- rename(Overlap, "IID"="f.eid")

SZ_merged <- merge(SZ_merged, Overlap, by = "IID", all.x = T)

SZ_merged <- SZ_merged %>% filter(is.na(COUNT))

SZ_score <- sapply(Scores_to_test, PES_PRS_SZ_UKBB)

## Extract beta, se, z, and p value for each gene

SZ_extract <- apply(SZ_score, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

SZ_results <- data.frame()

for (i in 1:length(SZ_extract)) {
  SZ_results <- rbind(SZ_results, SZ_extract[[i]])
}
rownames(SZ_results) <- Scores_to_test

write.csv(SZ_results, file="Results/Liberal_BIP_overlap_removed_SZ_UKBB_PES_PRS_results.csv", quote = F, row.names = T)



