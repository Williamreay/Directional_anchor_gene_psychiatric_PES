###############################

## Residualised PES - BIP TRAIN UKBB

## Penalised regression

##############################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_BIP_UKBB_TRAIN/")

library(dplyr)
library(data.table)
library(quantreg)
library(easyGgplot2)

set.seed(57373)

All_BIP_scores <- fread("All_BIP_TRAIN_UKBB_lassosum.txt", header = T, sep = " ")

Colnames_all_BIP_scores <- make.unique(names(All_BIP_scores))

colnames(All_BIP_scores) <- Colnames_all_BIP_scores

## Read in phenotype and covariate data

BIP_pheno <- fread("../../UKBB_phenotypes/BIP_phenotypes/BIP_TRAINING_UKBB_pheno.txt", header = T)

BIP_cov <- fread("../../UKBB_phenotypes/BIP_phenotypes/Covariates_TRAIN_BIP_UKBB.txt", header = T)

## Merge all three df

BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_BIP_scores, BIP_pheno, BIP_cov))

BIP_merged$Phenotype <- ifelse(BIP_merged$BIP == 1, "BIP", "HC")

BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)] <- lapply(BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)], function(x) c(scale(x)))

BIP_merged$Batch <- as.factor(BIP_merged$Batch)


## Kernel density plot for raw scores ##

ggplot2.density(data=BIP_merged, 
                xName='GRIN2A_PES', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of GRIN2A network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c("lightgreen",'#999999'))


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
                densityAlpha=0.3, densityLineSize=4, groupColors=c('lightgreen','#999999'))



Decile_resid_df <- BIP_merged %>%
  mutate_at(vars(ends_with("residuals")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("PES")), .funs = list(~ntile(.,10)))

table(Decile_resid_df$GRIN2A_residuals, Decile_resid_df$Phenotype)

table(Decile_resid_df$GRIN2A_PES, Decile_resid_df$Phenotype)

RPS17_decile_discrepancy <- Decile_resid_df %>% filter(RPS17_0_05 == 10 & RPS17_residuals < 10) %>%
  select(Phenotype, RPS17_0_05, RPS17_residuals) %>% filter(Phenotype == "BIP")


## Plot raw FES PES vs residualised FES PES

quantile(BIP_merged$GRIN2A_PES, probs = seq(0, 1, 0.1))

quantile(BIP_merged$GRIN2A_residuals, probs = seq(0, 1, 0.1))


BIP_merged$Plot_var_col <- as.factor(ifelse(BIP_merged$GRIN2A_PES > 1.30132773 & BIP_merged$GRIN2A_residuals > 1.29624489 , 1, 0))

ggplot(data = BIP_merged,
       aes(x = GRIN2A_PES, y=GRIN2A_residuals, colour = Plot_var_col)) +
  geom_point(alpha = 0.9, size = 1.4) +
  geom_vline(xintercept = 1.30132773 , linetype="longdash") +
  geom_hline(yintercept=1.29624489, linetype="longdash") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Raw GRIN2A network PES") +
  ylab("Residualised GRIN2A network PES") +
  scale_color_manual(name="Plot_var_col", values=c('#999999','turquoise3'))
