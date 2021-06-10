#############################

## Correlation of PES/PRS with self-reported MHQ phenotypes

## SZ and BIP excluded + controls in the discovery analyses

## Sensitivity analyses

## William Reay (2021)

############################

library(dplyr)
library(data.table)


setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_self_reported_mental_illness/")

## Read in PES/PRS

PES_PRS <- fread("../UKBB_MHQ_no_self_reported_biochem/Entire_EUR_genotyped_UKBB_scores/FINAL_SCORES_PES_PRS/All_SZ_BIP_PES_PRS_all_UKBB.txt", header = T)

PES_PRS <- PES_PRS[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28)]

## Load in MHQ self-report data

MHQ_self_report <- fread("UKBB_self_report_mental_illness_df_SZ_BIP_excl.txt", header = T)

MHQ_self_report <- rename(MHQ_self_report, "IID"="f.eid")

## Merge

Merged_MHQ_scoring <- merge(PES_PRS, MHQ_self_report, by = "IID")

Merged_MHQ_scoring$Batch <- as.factor(Merged_MHQ_scoring$Batch)

## Define list of scores

Scores_to_test <- as.list(colnames(Merged_MHQ_scoring[,2:15]))

## SZ CACNA1C PES and depression

Dep <- Merged_MHQ_scoring %>% filter(No_mental_health == 0 | Depression == 1)

Dep$SZ_CACNA1C_PES <- as.numeric(scale(Dep$SZ_CACNA1C_PES))
Dep$SZ_PRS <- as.numeric(scale(Dep$SZ_PRS))

Dep_test <- glm(Depression ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
  PC8 + PC9 + PC10 + Batch + SZ_CACNA1C_PES + SZ_PRS, family = "binomial",
  data = Dep)

anova(Dep_test, test="Chisq")

## SZ RPS17 PES and OCD

OCD <- Merged_MHQ_scoring %>% filter(No_mental_health == 0 | OCD == 1)

OCD$SZ_RSP17_PES <- as.numeric(scale(OCD$SZ_RPS17_PES))
OCD$SZ_PRS <- as.numeric(scale(OCD$SZ_PRS))

OCD_test <- glm(OCD ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                  PC8 + PC9 + PC10 + Batch + SZ_RPS17_PES + SZ_PRS, family = "binomial",
                data = OCD)

anova(OCD_test, test="Chisq")

## BIP GRIN2A PES and anxiety

Anx <- Merged_MHQ_scoring %>% filter(No_mental_health == 0 | Anx_nerves_or_GAD == 1)

Anx$BIP_GRIN2A_PES <- as.numeric(scale(Anx$BIP_GRIN2A_PES))
Anx$BIP_PRS <- as.numeric(scale(Anx$BIP_PRS))

BIP_test <- glm(Anx_nerves_or_GAD ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                  PC8 + PC9 + PC10 + Batch + BIP_GRIN2A_PES + BIP_PRS, family = "binomial",
                data = OCD)

anova(BIP_test, test="Chisq")
