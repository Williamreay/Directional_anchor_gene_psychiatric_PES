######################################

## BIP UKBB - PES/PRS clustering and thresholding

## William Reay (2021)

#########################################

set.seed(8764386)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/BIP_TRAIN_UKBB/")

library(mclust)
library(dplyr)
library(data.table)
library(factoextra)

All_BIP_scores <- fread("Combined_BIP_PES_PRS_train.txt", header = T, sep = " ")

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

GMM_df <- BIP_merged %>% filter(BIP == 1) %>% select(PRS_0_005, FES_0_005, GRIN2A_0_005, RPS17_0_05, CACNA1C_0_005,
                                                   PCCB_0_05, FADS1_1)

GMM_df[ ,c(1:7)] <- lapply(GMM_df[ ,c(1:7)], function(x) c(scale(x)))

BIP_UKBB_BIC <- mclustBIC(GMM_df, legendArgs=list(x="top", ncol=3))

plot(BIP_UKBB_BIC)

BIP_UKBB_mod <- Mclust(GMM_df, x = BIP_UKBB_BIC)

plot(BIP_UKBB_mod, what = "classification", dimens = c(1,2,3))

plot(BIP_UKBB_mod, what = "classification", dimens = c(1, 4, 5))

plot(BIP_UKBB_mod, what = "classification")

fviz_cluster(BIP_UKBB_mod, GMM_df,
             palette = c("grey", "#2E9FDF"), 
             geom = "point",
             ellipse.type = "norm", 
             ggtheme = theme_bw(),
             main="BIP PES/PRS [finite GMM, VVI, components = 2]",
             show.clust.cent = F
)


## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

BIP_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)] <- lapply(BIP_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)], function(x) c(scale(x)))

## Identify number of BIP cases with top decile PES

Decile_raw_df <- BIP_merged %>% select(IID, Sex, Age, BIP, PRS_0_005, FES_0_005, GRIN2A_0_005, RPS17_0_05, CACNA1C_0_005,
                                       PCCB_0_05, FADS1_1)

Decile_df <- Decile_raw_df %>%
  mutate_at(vars(ends_with("5")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("1")), .funs = list(~ntile(.,10)))

## Assign individuals in the top PES decile

Top_decile_PES <- Decile_df %>%
  mutate_at(vars(ends_with("5")), .funs = list(~if_else(. == 10,1, 0))) %>% mutate_at(vars(ends_with("1")), .funs = list(~if_else(. == 10,1, 0)))


Top_decile_PES$PES_decile_sum <- rowSums(Top_decile_PES[, c(6:11)])

Top_decile_PES$Binary_decile <- ifelse(Top_decile_PES$PES_decile_sum > 0, 1, 0)

table(Top_decile_PES$BIP, Top_decile_PES$PES_decile_sum)

Decile_BIP <- glm(BIP ~  Binary_decile,
                 family = "binomial", data = Top_decile_PES)

PRS_decile <- glm(BIP ~ PRS_0_005 + Binary_decile,
                  family = "binomial", data = Top_decile_PES)

## Identify individuals in the bottom decile of PRS

Bottom_decile_PRS <- Decile_raw_df %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))

Bottom_decile_PRS_only <- Bottom_decile_PRS %>% filter(PRS_0_005 == 1)

Bottom_top_decile_merged <- merge(Bottom_decile_PRS_only, Top_decile_PES, by = "IID")

table(Bottom_top_decile_merged$BIP.y, Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)

BIP_dec <- glm(BIP.x ~ Sex.x + Age.x + PES_decile_sum, family = "binomial", data = Bottom_top_decile_merged)

## Repeat above but only for BIP cases

BIP_Decile_df <- Decile_raw_df %>% filter(BIP == 1) %>%
  mutate_at(vars(ends_with("5")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("1")), .funs = list(~ntile(.,10)))

BIP_top_decile_PES <- BIP_Decile_df %>%
  mutate_at(vars(ends_with("5")), .funs = list(~if_else(. == 10,1, 0))) %>% mutate_at(vars(ends_with("1")), .funs = list(~if_else(. == 10,1, 0)))


BIP_top_decile_PES$PES_decile_sum <- rowSums(BIP_top_decile_PES[, c(6:11)])

BIP_top_decile_PES$Binary_decile <- ifelse(BIP_top_decile_PES$PES_decile_sum > 0, 1, 0)

## Identify individuals in the bottom decile of PRS

BIP_Bottom_decile_PRS <- Decile_raw_df %>% filter(BIP == 1) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))

BIP_Bottom_decile_PRS_only <- BIP_Bottom_decile_PRS %>% filter(PRS_0_005 == 1)

BIP_Bottom_top_decile_merged <- merge(BIP_Bottom_decile_PRS_only, BIP_top_decile_PES, by = "IID")

table(BIP_Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)
