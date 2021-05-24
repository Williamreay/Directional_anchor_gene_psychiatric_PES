#############################################

## Association of PES and PRS with BIP in the UKBB cohort - TRAIN

## Penalised regression

## BIP cases and double the number of controls (randomly selected from MHQ participants)

## William Reay (2021)

#############################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_BIP_UKBB_TRAIN/")

library(data.table)
library(dplyr)
library(rcompanion)
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(viridis)


## Read in concatenated PES/PRS file

All_BIP_scores <- fread("All_BIP_TRAIN_UKBB_lassosum.txt", header = T, sep = " ")

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

write.csv(BIP_results, file="Results/Lassosum_BIP_UKBB_PES_PRS_results.csv", quote = F, row.names = T)


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

write.csv(R2_liability_results, file="Results/Lassosum_BIP_Liability_scale_r2_PES_PRS.csv", row.names = T, quote = F)


## Test covariation of RPS17 PES with PRS at same threshold

RPS17_full <- glm(BIP ~ Sex + Age + PC1 + PC2 + PC3 +
                    PC4 + PC5 + Batch + PRS + RPS17_PES, family = "binomial", data = BIP_merged)

anova(RPS17_full, test = "Chisq")

RPS17_PRS_cor <- lm(PRS_0_05 ~ Sex + Age + PC1 + PC2 + PC3 +
                      PC4 + PC5 + Batch + RPS17_0_05, data = BIP_merged)

## Test covariation of GRIN2A PES with PRS at same threshold

GRIN2A_full <- glm(BIP ~ Sex + Age + PC1 + PC2 + PC3 +
                     PC4 + PC5 + Batch + PRS + GRIN2A_PES, family = "binomial", data = BIP_merged)

anova(GRIN2A_full, test = "Chisq")

GRIN2A_PRS_cor <- lm(PRS_0_005 ~ Sex + Age + PC1 + PC2 + PC3 +
                       PC4 + PC5 + Batch + GRIN2A_0_005, data = BIP_merged)


## Test covariation of PCCB PES with PRS at same threshold

PCCB_full <- glm(BIP ~ Sex + Age + PC1 + PC2 + PC3 +
                   PC4 + PC5 + Batch + PRS + PCCB_PES, family = "binomial", data = BIP_merged)

anova(PCCB_full, test = "Chisq")

PCCB_PRS_cor <- lm(PRS_0_05 ~ Sex + Age + PC1 + PC2 + PC3 +
                     PC4 + PC5 + Batch + PCCB_0_05, data = BIP_merged)

## Test covariation of CACNA1C PES with PRS at same threshold

CACNA1C_full <- glm(BIP ~ Sex + Age + PC1 + PC2 + PC3 +
                      PC4 + PC5 + Batch + PRS + CACNA1C_PES, family = "binomial", data = BIP_merged)

anova(CACNA1C_full, test = "Chisq")



## Correlation plot

Cor_BIP <- cor(BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)])

write.table(Cor_BIP, file="Results/BIP_train_corr_mat.txt", quote = F, row.names = T, col.names = T)

pval <- psych::corr.test(BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)], adjust = "none")$p

Sig_corr_plot <- corrplot(Cor_BIP, tl.cex = 0.5, p.mat = pval, insig = "blank",
                          sig.level = 0.05, type="upper", method="number",
                          tl.col="black", tl.srt=45)

Bon_corr_plot <- corrplot(Cor_BIP, tl.cex = 0.5, p.mat = pval, insig = "blank",
                          sig.level = 0.001020408, type="upper", method="number",
                          tl.col="black", tl.srt=45)

## Make forest plot of BIP OR per SD

BIP_results$Scores <- Scores_to_test


Plot_df <-  BIP_results %>% filter(Scores == "FES_PES" |
                                     Scores == "RPS17_PES" | Scores == "GRIN2A_PES" |
                                     Scores == "CACNA1C_PES" | Scores == "PCCB_PES" |
                                     Scores == "FADS1_PES")

Plot_df$Score_name <- c("CACNA1C network PES", "FADS1 network PES", 
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

FP_BIP + scale_color_viridis(discrete=TRUE)

## Remove cases who also meet a BIP criteria

Overlap <- fread("../../UKBB_phenotypes/IDs_meeting_criteria_for_SZ_and_BIP.txt", header = T)

Overlap <- rename(Overlap, "IID"="f.eid")

BIP_merged <- merge(BIP_merged, Overlap, by = "IID", all.x = T)

BIP_merged <- BIP_merged %>% filter(is.na(COUNT))

BIP_score <- sapply(Scores_to_test, PES_PRS_BIP_UKBB)

## Extract beta, se, z, and p value for each gene

BIP_extract <- apply(BIP_score, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

BIP_results <- data.frame()

for (i in 1:length(BIP_extract)) {
  BIP_results <- rbind(BIP_results, BIP_extract[[i]])
}
rownames(BIP_results) <- Scores_to_test

write.csv(BIP_results, file="Results/SZ_overlap_removed_BIP_UKBB_PES_PRS_results.csv", quote = F, row.names = T)
