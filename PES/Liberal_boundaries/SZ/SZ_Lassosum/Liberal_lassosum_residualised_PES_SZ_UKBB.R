###############################

## Residualised PES - SZ UKBB

##############################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/")

library(dplyr)
library(data.table)
library(quantreg)
library(easyGgplot2)

set.seed(57373)

All_SZ_scores <- fread("All_SZ_PES_PRS_validation_UKBB_PES_PRS.txt", header = T, sep = " ")

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores

## Read in phenotype and covariate data

SZ_pheno <- fread("../../UKBB_phenotypes/SZ_phenotypes/SZ_UKBB_pheno.txt", header = T)

SZ_cov <- fread("../../UKBB_phenotypes/SZ_phenotypes/Covariates_SZ_UKBB.txt", header = T)

## Merge all three df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_SZ_scores, SZ_pheno, SZ_cov))

SZ_merged$Phenotype <- ifelse(SZ_merged$SZ == 1, "SZ", "HC")

SZ_merged[,c(5, 10, 15, 20, 25, 30, 35)] <- lapply(SZ_merged[,c(5, 10, 15, 20, 25, 30, 35)], function(x) c(scale(x)))

SZ_merged$Batch <- as.factor(SZ_merged$Batch)


## Kernel density plot for raw scores ##

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

## Residualise the FES and GRIN2A PES

ggplot2.density(data=SZ_merged, 
                xName='GRIN2A_residuals', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of residualised GRIN2A network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c('#999999','dodgerblue1'))


## Plot raw GRIN2A PES vs residualised GRIN2A PES

quantile(SZ_merged$GRIN2A_PES, probs = seq(0, 1, 0.1))

quantile(SZ_merged$GRIN2A_residuals, probs = seq(0, 1, 0.1))


SZ_merged$Plot_var_col <- as.factor(ifelse(SZ_merged$GRIN2A_PES > 1.26710161 & SZ_merged$GRIN2A_residuals > 1.26504500, 1, 0))

SZ_merged$Plot_var_col_PRS <- as.factor(ifelse(SZ_merged$GRIN2A_PES > 1.26710161 & SZ_merged$PRS < -1.27575454, 1, 0))

ggplot(data = SZ_merged,
       aes(x = GRIN2A_PES, y=GRIN2A_residuals, colour = Plot_var_col)) +
  geom_point(alpha = 0.9, size = 1.4) +
  geom_vline(xintercept = 1.26504500 , linetype="longdash") +
  geom_hline(yintercept=1.26710161, linetype="longdash") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Raw GRIN2A network PES") +
  ylab("Residualised GRIN2A network PES") +
  scale_color_manual(name="Plot_var_col", values=c('#999999','slateblue1'))

ggplot(data = SZ_merged,
       aes(x = GRIN2A_PES, y=PRS, colour = Plot_var_col_PRS)) +
  geom_point(alpha = 0.9, size = 1.4) +
  geom_vline(xintercept = 1.26504500 , linetype="longdash") +
  geom_hline(yintercept=-1.27575454, linetype="longdash") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Raw GRIN2A network PES") +
  ylab("Genome-wide PRS") +
  scale_color_manual(name="Plot_var_col", values=c('#999999','sienna1'))
