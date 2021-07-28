######################################

## BIP UKBB - PES/PRS clustering and thresholding + residualised PES

## William Reay (2021)

#########################################

set.seed(963632)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/")

library(mclust)
library(dplyr)
library(data.table)
library(easyGgplot2)
library(corrplot)
library(vcd)

All_BIP_scores <- fread("Best_DA_gene_PES/BIP_best_con_or_lib/BIP_UKBB_TRAIN_score_files/Combined_BIP_UKBB_TRAIN_best_scores.txt", header = T, sep = " ")

Colnames_all_BIP_scores <- make.unique(names(All_BIP_scores))

colnames(All_BIP_scores) <- Colnames_all_BIP_scores

## Read in phenotype and covariate data

BIP_pheno <- fread("UKBB_phenotypes/BIP_phenotypes/BIP_TRAINING_UKBB_pheno.txt", header = T)

BIP_cov <- fread("UKBB_phenotypes/BIP_phenotypes/Covariates_TRAIN_BIP_UKBB.txt", header = T)

## Merge all three df

BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_BIP_scores, BIP_pheno, BIP_cov))

BIP_merged$Batch <- as.factor(BIP_merged$Batch)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

BIP_merged[,c(2, 4, 6, 8, 10, 12, 14)] <- lapply(BIP_merged[,c(2, 4, 6, 8, 10, 12, 14)], function(x) c(scale(x)))

## Define function to adjust each PES for PRS - for significant PES - PRS pairs perform a chi-square test of resid dev

Scores_to_test <- as.list(colnames(BIP_merged[,c(2, 4, 6, 8, 10, 14)]))

PRS_cov <- function(v){
  Score_model <- glm(glue::glue('BIP ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + {v} + Batch + PRS'), family = "binomial", data = BIP_merged)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

BIP_score <- sapply(Scores_to_test, PRS_cov)

## Extract beta, se, z, and p value for each gene

BIP_extract <- apply(BIP_score, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

BIP_results <- data.frame()

for (i in 1:length(BIP_extract)) {
  BIP_results <- rbind(BIP_results, BIP_extract[[i]])
}
rownames(BIP_results) <- Scores_to_test

write.csv(BIP_results, file="Best_DA_gene_PES/BIP_best_con_or_lib/Covariation_for_PRS_UKBB_thresholding_BIP/Covaried_for_PRS_BIP_PES.csv", quote = F, row.names = T)

BIP_merged <- rename(BIP_merged, "CACNA1C PES"="CACNA1C_PES",
                    "FADS1 PES" = "FADS1_1", "GRIN2A PES"="GRIN2A_PES",
                    "PCCB PES"="PCCB_PES", "RPS17 PES"="RPS17_PES", "FES PES"="FES_PES")

## Construct correlation plot

Cor_BIP <- cor(BIP_merged[,c(2, 4, 6, 8, 10, 12, 14)])

pval <- psych::corr.test(BIP_merged[,c(2, 4, 6, 8, 10, 12, 14)], adjust = "none")$p

Sig_corr_plot <- corrplot(Cor_BIP, tl.cex = 0.5, p.mat = pval, insig = "blank",
                          sig.level = 0.05, type="upper", method="number",
                          tl.col="black", tl.srt=45)

## Identify number of cases with top decile PES

Decile_df <- BIP_merged %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("PRS")), .funs = list(~ntile(.,10)))

## Assign individuals in the top PES decile

Top_decile_PES <- Decile_df %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~if_else(. == 10,1, 0)))

Top_decile_PES$PES_decile_sum <- rowSums(Top_decile_PES[, c(2, 4, 6, 8, 10, 14)])

Top_decile_PES$Binary_decile <- ifelse(Top_decile_PES$PES_decile_sum > 0, 1, 0)

table(Top_decile_PES$BIP, Top_decile_PES$PES_decile_sum)

Dec_tab <- table(Top_decile_PES$BIP, Top_decile_PES$Binary_decile)

Decile_BIP <- glm(BIP ~  Binary_decile + PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                   PC8 + PC9 + PC10 + Batch + Sex + Age,
                 family = "binomial", data = Top_decile_PES)

PRS_decile <- glm(Binary_decile ~ PRS,
                  family = "binomial", data = Top_decile_PES)

## Make mosaic plot of BIP cases and controls per decile grouping

mosaic( ~ BIP + Binary_decile, data = Top_decile_PES,
        highlighting = "BIP", highlighting_fill = c("grey", "plum"))


## Identify individuals in the bottom decile of PRS

Bottom_decile_PRS <- Decile_df %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))

Bottom_decile_PRS_only <- Bottom_decile_PRS %>% filter(PRS == 1)

Bottom_top_decile_merged <- merge(Bottom_decile_PRS_only, Top_decile_PES, by = "IID")

table(Bottom_top_decile_merged$BIP.y, Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)

BIP_dec <- glm(BIP.x ~ Sex.x + Age.x + Binary_decile, family = "binomial", data = Bottom_top_decile_merged)

## Repeat above but only for BIP cases

BIP_Decile_df <- Decile_df %>% filter(BIP == 1) %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("PRS")), .funs = list(~ntile(.,10)))

BIP_top_decile_PES <- BIP_Decile_df %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~if_else(. == 10,1, 0)))


BIP_top_decile_PES$PES_decile_sum <- rowSums(BIP_top_decile_PES[, c(2, 4, 6, 8, 10, 14)])

BIP_top_decile_PES$Binary_decile <- ifelse(BIP_top_decile_PES$PES_decile_sum > 0, 1, 0)

## Identify individuals in the bottom decile of PRS

BIP_Bottom_decile_PRS <- Decile_df %>% filter(BIP == 1) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))

BIP_Bottom_decile_PRS_only <- BIP_Bottom_decile_PRS %>% filter(PRS == 1)

BIP_Bottom_top_decile_merged <- merge(BIP_Bottom_decile_PRS_only, BIP_top_decile_PES, by = "IID")

table(BIP_Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)

## Residualised GRIN2A PES

BIP_merged$Phenotype <- ifelse(BIP_merged$BIP == 1, "BIP", "HC")

ggplot2.density(data=BIP_merged, 
                xName='GRIN2A_PES', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of GRIN2A network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c('#999999','slateblue1'))

## Residualised PES ##

GRIN2A_resid_mod <- lm(GRIN2A_PES ~ PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS,
                       data = BIP_merged)

GRIN2A_residuals <- as.data.frame(GRIN2A_resid_mod$residuals)

colnames(GRIN2A_residuals)[1] <- "GRIN2A_residuals"

BIP_merged$GRIN2A_residuals <- as.data.frame(scale(GRIN2A_residuals$GRIN2A_residuals))

ggplot2.density(data=BIP_merged, 
                xName='GRIN2A_residuals', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of residualised GRIN2A network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c('#999999','slateblue1'))


Decile_resid_df <- BIP_merged %>%
  mutate_at(vars(ends_with("residuals")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~ntile(.,10)))

table(Decile_resid_df$GRIN2A_residuals, Decile_resid_df$Phenotype)

table(Decile_resid_df$GRIN2A_PES, Decile_resid_df$Phenotype)

GRIN2A_decile_discrepancy <- Decile_resid_df %>% filter(GRIN2A_PES == 10 & GRIN2A_residuals < 9)

table(Decile_resid_df$GRIN2A_PES, Decile_resid_df$Phenotype)

table(Decile_resid_df$GRIN2A_residuals, Decile_resid_df$Phenotype)

## Plot raw FES PES vs residualised FES PES

quantile(BIP_merged$GRIN2A_PES, probs = seq(0, 1, 0.1))

quantile(BIP_merged$GRIN2A_residuals, probs = seq(0, 1, 0.1))

quantile(BIP_merged$PRS, probs = seq(0, 1, 0.1))

BIP_merged$Plot_var_col <- as.factor(ifelse(BIP_merged$GRIN2A_PES > 1.31265716 & BIP_merged$GRIN2A_residuals > 1.29474340, 1, 0))

BIP_merged$Plot_var_col_2 <- as.factor(ifelse(BIP_merged$GRIN2A_PES > 1.31265716 & BIP_merged$PRS < -1.27051916 , 1, 0))

ggplot(data = BIP_merged,
       aes(x = GRIN2A_PES, y=GRIN2A_residuals, colour = Plot_var_col)) +
  geom_point(alpha = 0.9, size = 1.4) +
  geom_vline(xintercept = 1.31039008 , linetype="longdash") +
  geom_hline(yintercept=1.29834122, linetype="longdash") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Raw GRIN2A network PES") +
  ylab("Residualised GRIN2A network PES") +
  scale_color_manual(name="Plot_var_col", values=c('#999999','coral2'))

## PRS vs PES

ggplot(data = BIP_merged,
       aes(x = GRIN2A_PES, y=PRS, colour = Plot_var_col_2)) +
  geom_point(alpha = 0.9, size = 1.4) +
  geom_vline(xintercept = 1.31039008 , linetype="longdash") +
  geom_hline(yintercept=-1.26673774, linetype="longdash") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Raw GRIN2A network PES") +
  ylab("Genome-wide PRS") +
  scale_color_manual(name="Plot_var_col_2", values=c('#999999','cornflowerblue'))
