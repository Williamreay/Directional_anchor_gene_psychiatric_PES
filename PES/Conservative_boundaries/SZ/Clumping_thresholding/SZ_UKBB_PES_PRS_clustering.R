######################################

## SZ UKBB - PES/PRS clustering and thresholding

## William Reay (2021)

#########################################

set.seed(8764386)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/UKBB_SZ/")

library(mclust)
library(dplyr)
library(data.table)

All_SZ_scores <- fread("Combined_SZ_all_PES_PRS_UKBB.txt", header = T, sep = " ")

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores

## Read in phenotype and covariate data

SZ_pheno <- fread("../../UKBB_phenotypes/SZ_phenotypes/SZ_UKBB_pheno.txt", header = T)

SZ_cov <- fread("../../UKBB_phenotypes/SZ_phenotypes/Covariates_SZ_UKBB.txt", header = T)

## Merge all three df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_SZ_scores, SZ_pheno, SZ_cov))

SZ_merged$Batch <- as.factor(SZ_merged$Batch)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

SZ_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)] <- lapply(SZ_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)], function(x) c(scale(x)))


## Select the most variance explained thresholds and SZ cases only

GMM_df <- SZ_merged %>% filter(SZ == 1) %>% select(PRS_0_05, FES_0_5, GRIN2A_0_005, RPS17_1, CACNA1C_0_005,
                               PCCB_0_05, FADS1_0_05)

SZ_UKBB_BIC <- mclustBIC(GMM_df, legendArgs=list(x="top", ncol=3))

plot(SZ_UKBB_BIC)

SZ_UKBB_mod <- Mclust(GMM_df, x = SZ_UKBB_BIC)

plot(SZ_UKBB_mod, what = "classification", dimens = c(1,2,3))

plot(SZ_UKBB_mod, what = "classification", dimens = c(1, 4, 5, 6))


## Identify number of SZ cases with top decile PES

Decile_raw_df <- SZ_merged %>% select(IID, Sex, Age, SZ, PRS_0_05, FES_0_5, GRIN2A_0_005, RPS17_1, CACNA1C_0_005,
                              PCCB_0_05, FADS1_0_05)

Decile_df <- Decile_raw_df %>%
  mutate_at(vars(ends_with("5")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("1")), .funs = list(~ntile(.,10)))

## Assign individuals in the top PES decile

Top_decile_PES <- Decile_df %>%
  mutate_at(vars(ends_with("5")), .funs = list(~if_else(. == 10,1, 0))) %>% mutate_at(vars(ends_with("1")), .funs = list(~if_else(. == 10,1, 0)))


Top_decile_PES$PES_decile_sum <- rowSums(Top_decile_PES[, c(6:11)])

Top_decile_PES$Binary_decile <- ifelse(Top_decile_PES$PES_decile_sum > 0, 1, 0)

table(Top_decile_PES$SZ, Top_decile_PES$PES_decile_sum)

Decile_SZ <- glm(SZ ~  Binary_decile + PRS_0_05,
                 family = "binomial", data = Top_decile_PES)

PRS_decile <- glm(Binary_decile ~ PRS_0_05,
                  family = "binomial", data = Top_decile_PES)

## Identify individuals in the bottom decile of PRS

Bottom_decile_PRS <- Decile_raw_df %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))
  
Bottom_decile_PRS_only <- Bottom_decile_PRS %>% filter(PRS_0_05 == 1)

Bottom_top_decile_merged <- merge(Bottom_decile_PRS_only, Top_decile_PES, by = "IID")

table(Bottom_top_decile_merged$SZ.y, Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)

SZ_dec <- glm(SZ.x ~ Sex.x + Age.x + PES_decile_sum, family = "binomial", data = Bottom_top_decile_merged)

## Repeat above but only for SZ cases

SZ_Decile_df <- Decile_raw_df %>% filter(SZ == 1) %>%
  mutate_at(vars(ends_with("5")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("1")), .funs = list(~ntile(.,10)))

SZ_top_decile_PES <- SZ_Decile_df %>%
  mutate_at(vars(ends_with("5")), .funs = list(~if_else(. == 10,1, 0))) %>% mutate_at(vars(ends_with("1")), .funs = list(~if_else(. == 10,1, 0)))


SZ_top_decile_PES$PES_decile_sum <- rowSums(SZ_top_decile_PES[, c(6:11)])

SZ_top_decile_PES$Binary_decile <- ifelse(SZ_top_decile_PES$PES_decile_sum > 0, 1, 0)

## Identify individuals in the bottom decile of PRS

SZ_Bottom_decile_PRS <- Decile_raw_df %>% filter(SZ == 1) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))

SZ_Bottom_decile_PRS_only <- SZ_Bottom_decile_PRS %>% filter(PRS_0_05 == 1)

SZ_Bottom_top_decile_merged <- merge(SZ_Bottom_decile_PRS_only, SZ_top_decile_PES, by = "IID")

table(SZ_Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)
