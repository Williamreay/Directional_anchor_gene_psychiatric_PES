###############################

## Residualised PES - SZ UKBB

##############################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/UKBB_SZ/")

library(dplyr)
library(data.table)
library(quantreg)
library(easyGgplot2)

set.seed(57373)

All_SZ_scores <- fread("Combined_SZ_all_PES_PRS_UKBB.txt", header = T, sep = " ")

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores

## Read in phenotype and covariate data

SZ_pheno <- fread("../../UKBB_phenotypes/SZ_phenotypes/SZ_UKBB_pheno.txt", header = T)

SZ_cov <- fread("../../UKBB_phenotypes/SZ_phenotypes/Covariates_SZ_UKBB.txt", header = T)

## Merge all three df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_SZ_scores, SZ_pheno, SZ_cov))

SZ_merged$Phenotype <- ifelse(SZ_merged$SZ == 1, "SZ", "HC")

SZ_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)] <- lapply(SZ_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)], function(x) c(scale(x)))

SZ_merged$Batch <- as.factor(SZ_merged$Batch)


## Kernel density plot for raw scores ##

ggplot2.density(data=SZ_merged, 
                xName='FES_0_5', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of FES network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c('#999999','#E69F00'))

ggplot2.density(data=SZ_merged, 
                xName='GRIN2A_0_005', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of GRIN2A network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c('#999999','dodgerblue1'))


## Residualised PES ##

FES_resid_mod <- lm(FES_0_5 ~ PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS_0_5,
                    data = SZ_merged)

GRIN2A_resid_mod <- lm(GRIN2A_0_005 ~ PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS_0_005,
                       data = SZ_merged)

FES_residuals <- as.data.frame(FES_resid_mod$residuals)

GRIN2A_residuals <- as.data.frame(GRIN2A_resid_mod$residuals)

colnames(FES_residuals)[1] <- "FES_residuals"
colnames(GRIN2A_residuals)[1] <- "GRIN2A_residuals"

SZ_merged$FES_residuals <- as.data.frame(scale(FES_residuals$FES_residuals))
SZ_merged$GRIN2A_residuals <- as.data.frame(scale(GRIN2A_residuals$GRIN2A_residuals))

## Residualise the FES and GRIN2A PES

ggplot2.density(data=SZ_merged, 
                xName='FES_residuals', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of residualised FES network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c('#999999','#E69F00'))

ggplot2.density(data=SZ_merged, 
                xName='GRIN2A_residuals', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of residualised GRIN2A network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c('#999999','dodgerblue1'))


Decile_resid_df <- SZ_merged %>%
  mutate_at(vars(ends_with("residuals")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("5")), .funs = list(~ntile(.,10)))

table(Decile_resid_df$FES_residuals, Decile_resid_df$Phenotype)

table(Decile_resid_df$FES_0_5, Decile_resid_df$Phenotype)

FES_decile_discrepancy <- Decile_resid_df %>% filter(FES_0_5 == 10 & FES_residuals < 9) %>%
  select(Phenotype, FES_0_5, FES_residuals)

table(Decile_resid_df$GRIN2A_0_005, Decile_resid_df$Phenotype)

table(Decile_resid_df$GRIN2A_residuals, Decile_resid_df$Phenotype)

## Plot raw FES PES vs residualised FES PES

quantile(SZ_merged$FES_0_5, probs = seq(0, 1, 0.1))

quantile(SZ_merged$FES_residuals, probs = seq(0, 1, 0.1))


SZ_merged$Plot_var_col <- as.factor(ifelse(SZ_merged$FES_0_5 > 1.26282143 & SZ_merged$FES_residuals > 1.27726849, 1, 0))

ggplot(data = SZ_merged,
       aes(x = FES_0_5, y=FES_residuals, colour = Plot_var_col)) +
  geom_point(alpha = 0.9, size = 1.4) +
  geom_vline(xintercept = 1.26282143 , linetype="longdash") +
  geom_hline(yintercept=1.27726849, linetype="longdash") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Raw FES network PES") +
  ylab("Residualised PES network PES") +
  scale_color_manual(name="Plot_var_col", values=c('#999999','dodgerblue1'))