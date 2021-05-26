######################################

## SZ UKBB - PES/PRS clustering and thresholding + residualised PES

## William Reay (2021)

#########################################

set.seed(8764386)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/")

library(mclust)
library(dplyr)
library(data.table)
library(easyGgplot2)
library(vcd)

All_SZ_scores <- fread("Best_DA_gene_PES/SZ_best_con_or_lib/SZ_UKBB_TRAIN_score_files/Best_SZ_PES_PRS_UKBB_combined.txt", header = T, sep = " ")

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores

## Read in phenotype and covariate data

SZ_pheno <- fread("UKBB_phenotypes/SZ_phenotypes/SZ_UKBB_pheno.txt", header = T)

SZ_cov <- fread("UKBB_phenotypes/SZ_phenotypes/Covariates_SZ_UKBB.txt", header = T)

## Merge all three df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_SZ_scores, SZ_pheno, SZ_cov))

SZ_merged$Batch <- as.factor(SZ_merged$Batch)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

SZ_merged[,c(2, 4, 6, 11, 16, 21, 26)] <- lapply(SZ_merged[,c(2, 4, 6, 11, 16, 21, 26)], function(x) c(scale(x)))

## Define function to adjust each PES for PRS - for significant PES - PRS pairs perform a chi-square test of resid dev

Scores_to_test <- as.list(colnames(SZ_merged[,c(2, 4, 6, 11, 16, 21)]))

PRS_cov <- function(v){
  Score_model <- glm(glue::glue('SZ ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + {v} + Batch + PRS'), family = "binomial", data = SZ_merged)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

SZ_score <- sapply(Scores_to_test, PRS_cov)

## Extract beta, se, z, and p value for each gene

SZ_extract <- apply(SZ_score, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

SZ_results <- data.frame()

for (i in 1:length(SZ_extract)) {
  SZ_results <- rbind(SZ_results, SZ_extract[[i]])
}
rownames(SZ_results) <- Scores_to_test

write.csv(SZ_results, file="Best_DA_gene_PES/SZ_best_con_or_lib/Covariation_for_PRS_UKBB_thresholding/Covaried_for_PRS_SZ_PES.csv", quote = F, row.names = T)

SZ_merged <- rename(SZ_merged, "CACNA1C PES"="CACNA1C_0_005",
                    "FADS1 PES" = "FADS1_0_05", "GRIN2A PES"="GRIN2A_PES",
                    "PCCB PES"="PCCB_PES", "RPS17 PES"="RPS17_PES", "FES PES"="FES_0_05")

## Construct correlation plot

Cor_SZ <- cor(SZ_merged[,c(2, 4, 6, 11, 16, 21, 26)])

pval <- psych::corr.test(SZ_merged[,c(2, 4, 6, 11, 16, 21, 26)], adjust = "none")$p

Sig_corr_plot <- corrplot(Cor_SZ, tl.cex = 0.5, p.mat = pval, insig = "blank",
                          sig.level = 0.05, type="upper", method="number",
                          tl.col="black", tl.srt=45)

## Identify number of cases with top decile PES

Decile_df <- SZ_merged %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("PRS")), .funs = list(~ntile(.,10)))

## Assign individuals in the top PES decile

Top_decile_PES <- Decile_df %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~if_else(. == 10,1, 0)))

Top_decile_PES$PES_decile_sum <- rowSums(Top_decile_PES[, c(2, 4, 6, 11, 16, 21)])

Top_decile_PES$Binary_decile <- ifelse(Top_decile_PES$PES_decile_sum > 0, 1, 0)

table(Top_decile_PES$SZ, Top_decile_PES$PES_decile_sum)

Decile_SZ <- glm(SZ ~  Binary_decile + PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                   PC8 + PC9 + PC10 + Batch + Sex + Age,
                 family = "binomial", data = Top_decile_PES)

PRS_decile <- glm(Binary_decile ~ PRS,
                  family = "binomial", data = Top_decile_PES)

## Make mosaic plot of BIP cases and controls per decile grouping

mosaic(~ SZ + Binary_decile, data = Top_decile_PES,
        highlighting = "SZ", highlighting_fill = c("grey", "dodgerblue"))

## Identify individuals in the bottom decile of PRS

Bottom_decile_PRS <- Decile_df %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))

Bottom_decile_PRS_only <- Bottom_decile_PRS %>% filter(PRS == 1)

Bottom_top_decile_merged <- merge(Bottom_decile_PRS_only, Top_decile_PES, by = "IID")

table(Bottom_top_decile_merged$SZ.y, Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)

SZ_dec <- glm(SZ.x ~ Sex.x + Age.x + Binary_decile, family = "binomial", data = Bottom_top_decile_merged)

## Repeat above but only for SZ cases

SZ_Decile_df <- Decile_df %>% filter(SZ == 1) %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("PRS")), .funs = list(~ntile(.,10)))

SZ_top_decile_PES <- SZ_Decile_df %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~if_else(. == 10,1, 0)))


SZ_top_decile_PES$PES_decile_sum <- rowSums(SZ_top_decile_PES[, c(2, 4, 6, 11, 16, 21)])

SZ_top_decile_PES$Binary_decile <- ifelse(SZ_top_decile_PES$PES_decile_sum > 0, 1, 0)

## Identify individuals in the bottom decile of PRS

SZ_Bottom_decile_PRS <- Decile_df %>% filter(SZ == 1) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))

SZ_Bottom_decile_PRS_only <- SZ_Bottom_decile_PRS %>% filter(PRS == 1)

SZ_Bottom_top_decile_merged <- merge(SZ_Bottom_decile_PRS_only, SZ_top_decile_PES, by = "IID")

table(SZ_Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)

## Residualised GRIN2A PES

SZ_merged$Phenotype <- ifelse(SZ_merged$SZ == 1, "SZ", "HC")

ggplot2.density(data=SZ_merged, 
                xName='GRIN2A_PES', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of GRIN2A network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c('#999999','dodgerblue1'))

## Residualised PES ##

GRIN2A_resid_mod <- lm(GRIN2A_PES ~ PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS,
                       data = SZ_merged)

GRIN2A_residuals <- as.data.frame(GRIN2A_resid_mod$residuals)

colnames(GRIN2A_residuals)[1] <- "GRIN2A_residuals"

SZ_merged$GRIN2A_residuals <- as.data.frame(scale(GRIN2A_residuals$GRIN2A_residuals))

ggplot2.density(data=SZ_merged, 
                xName='GRIN2A_residuals', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of residualised GRIN2A network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c('#999999','dodgerblue1'))


Decile_resid_df <- SZ_merged %>%
  mutate_at(vars(ends_with("residuals")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~ntile(.,10)))

table(Decile_resid_df$GRIN2A_residuals, Decile_resid_df$Phenotype)

table(Decile_resid_df$GRIN2A_PES, Decile_resid_df$Phenotype)

GRIN2A_decile_discrepancy <- Decile_resid_df %>% filter(GRIN2A_PES == 10 & GRIN2A_residuals < 9)

table(Decile_resid_df$GRIN2A_PES, Decile_resid_df$Phenotype)

table(Decile_resid_df$GRIN2A_residuals, Decile_resid_df$Phenotype)

## Plot raw FES PES vs residualised FES PES

quantile(SZ_merged$GRIN2A_PES, probs = seq(0, 1, 0.1))

quantile(SZ_merged$GRIN2A_residuals, probs = seq(0, 1, 0.1))

quantile(SZ_merged$PRS, probs = seq(0, 1, 0.1))

SZ_merged$Plot_var_col <- as.factor(ifelse(SZ_merged$GRIN2A_PES > 1.27672190 & SZ_merged$GRIN2A_residuals > 1.248529572, 1, 0))

SZ_merged$Plot_var_col_2 <- as.factor(ifelse(SZ_merged$GRIN2A_PES > 1.27672190 & SZ_merged$PRS < -1.27575454, 1, 0))

ggplot(data = SZ_merged,
       aes(x = GRIN2A_PES, y=GRIN2A_residuals, colour = Plot_var_col)) +
  geom_point(alpha = 0.9, size = 1.4) +
  geom_vline(xintercept = 1.27672190 , linetype="longdash") +
  geom_hline(yintercept=1.248529572, linetype="longdash") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Raw GRIN2A network PES") +
  ylab("Residualised GRIN2A network PES") +
  scale_color_manual(name="Plot_var_col", values=c('#999999','#E69F00'))

## PRS vs PES

ggplot(data = SZ_merged,
       aes(x = GRIN2A_PES, y=PRS, colour = Plot_var_col_2)) +
  geom_point(alpha = 0.9, size = 1.4) +
  geom_vline(xintercept = 1.27672190 , linetype="longdash") +
  geom_hline(yintercept=-1.27575454, linetype="longdash") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Raw GRIN2A network PES") +
  ylab("Genome-wide PRS") +
  scale_color_manual(name="Plot_var_col_2", values=c('#999999','mediumpurple1'))

