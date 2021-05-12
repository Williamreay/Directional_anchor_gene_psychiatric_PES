######################################

## BIP UKBB - PES/PRS clustering and thresholding

## Penalised regression

## William Reay (2021)

#########################################

set.seed(78484)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_BIP_UKBB_TRAIN/")

library(mclust)
library(dplyr)
library(data.table)
library(factoextra)
library(ggplot2)
library(reshape2)

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




## Select the most variance explained thresholds

GMM_df <- BIP_merged %>% filter(BIP == 1) %>% select(PRS, FES_PES, GRIN2A_PES, RPS17_PES, CACNA1C_PES,
                                                     PCCB_PES, FADS1_PES)

GMM_df[ ,c(1:7)] <- lapply(GMM_df[ ,c(1:7)], function(x) c(scale(x)))

BIP_UKBB_BIC <- mclustBIC(GMM_df, legendArgs=list(x="top", ncol=3))

plot(BIP_UKBB_BIC)

BIP_UKBB_mod <- Mclust(GMM_df, x = BIP_UKBB_BIC)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)] <- lapply(BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)], function(x) c(scale(x)))

## Identify number of BIP cases with top decile PES

Decile_raw_df <- BIP_merged %>% select(IID, Sex, Age, BIP, PRS, FES_PES, GRIN2A_PES, RPS17_PES, CACNA1C_PES,
                                       PCCB_PES, FADS1_PES)

Decile_df <- Decile_raw_df %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("PRS")), .funs = list(~ntile(.,10)))

## Assign individuals in the top PES decile

Top_decile_PES <- Decile_df %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~if_else(. == 10,1, 0))) %>% mutate_at(vars(ends_with("1")), .funs = list(~if_else(. == 10,1, 0)))


Top_decile_PES$PES_decile_sum <- rowSums(Top_decile_PES[, c(6:11)])

Top_decile_PES$Binary_decile <- ifelse(Top_decile_PES$PES_decile_sum > 0, 1, 0)

table(Top_decile_PES$BIP, Top_decile_PES$PES_decile_sum)

Decile_BIP <- glm(BIP ~  Binary_decile,
                  family = "binomial", data = Top_decile_PES)

PRS_decile <- glm(BIP ~ PRS + Binary_decile,
                  family = "binomial", data = Top_decile_PES)

## Identify individuals in the bottom decile of PRS

Bottom_decile_PRS <- Decile_raw_df %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))

Bottom_decile_PRS_only <- Bottom_decile_PRS %>% filter(PRS == 1)

Bottom_top_decile_merged <- merge(Bottom_decile_PRS_only, Top_decile_PES, by = "IID")

table(Bottom_top_decile_merged$BIP.y, Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)

BIP_dec <- glm(BIP.x ~ Sex.x + Age.x + PES_decile_sum, family = "binomial", data = Bottom_top_decile_merged)

## Repeat above but only for BIP cases

BIP_Decile_df <- Decile_raw_df %>% filter(BIP == 1) %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("PRS")), .funs = list(~ntile(.,10)))

BIP_top_decile_PES <- BIP_Decile_df %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~if_else(. == 10,1, 0))) %>% mutate_at(vars(ends_with("1")), .funs = list(~if_else(. == 10,1, 0)))


BIP_top_decile_PES$PES_decile_sum <- rowSums(BIP_top_decile_PES[, c(6:11)])

BIP_top_decile_PES$Binary_decile <- ifelse(BIP_top_decile_PES$PES_decile_sum > 0, 1, 0)

## Identify individuals in the bottom decile of PRS

BIP_Bottom_decile_PRS <- Decile_raw_df %>% filter(BIP == 1) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))

BIP_Bottom_decile_PRS_only <- BIP_Bottom_decile_PRS %>% filter(PRS == 1)

BIP_Bottom_top_decile_merged <- merge(BIP_Bottom_decile_PRS_only, BIP_top_decile_PES, by = "IID")

table(BIP_Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)

## Make boxplot

Boxplot <- melt(BIP_merged,measure.vars=c('GRIN2A_PES','FES_PES','CACNA1C_PES', 'RPS17_PES', 'PCCB_PES', 'FADS1_PES', 'PRS'), id.vars='BIP')

Boxplot$BIP <- ifelse(Boxplot$BIP == 1, "BIP", "HC")

P <- ggplot(Boxplot) +
  geom_boxplot(aes(x=BIP, y=value, fill=variable)) +
  theme_bw() +
  xlab("Phenotype") +
  ylab("Score") +
  labs(fill = "Score name") +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype = "longdash")
  
