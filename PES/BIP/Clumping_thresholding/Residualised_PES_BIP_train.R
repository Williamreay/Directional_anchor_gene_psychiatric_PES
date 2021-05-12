###############################

## Residualised PES - BIP TRAIN UKBB

##############################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/BIP_TRAIN_UKBB/")

library(dplyr)
library(data.table)
library(quantreg)
library(easyGgplot2)

set.seed(57373)

All_BIP_scores <- fread("Combined_BIP_PES_PRS_train.txt", header = T, sep = " ")

Colnames_all_BIP_scores <- make.unique(names(All_BIP_scores))

colnames(All_BIP_scores) <- Colnames_all_BIP_scores

## Read in phenotype and covariate data

BIP_pheno <- fread("../../UKBB_phenotypes/BIP_phenotypes/BIP_TRAINING_UKBB_pheno.txt", header = T)

BIP_cov <- fread("../../UKBB_phenotypes/BIP_phenotypes/Covariates_TRAIN_BIP_UKBB.txt", header = T)

## Merge all three df

BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_BIP_scores, BIP_pheno, BIP_cov))

BIP_merged$Phenotype <- ifelse(BIP_merged$BIP == 1, "BIP", "HC")

BIP_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)] <- lapply(BIP_merged[,c(3:6, 9:12, 15:18, 21:24, 27:30, 33:36, 39:42)], function(x) c(scale(x)))

BIP_merged$Batch <- as.factor(BIP_merged$Batch)


## Kernel density plot for raw scores ##

ggplot2.density(data=BIP_merged, 
                xName='RPS17_0_05', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of RPS17 network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c("darkorchid1",'#999999'))


## Residualised PES ##

RPS17_resid_mod <- lm(RPS17_0_05 ~ PC1 + PC2 + PC3 + PC4 + PC5 + Batch + PRS_0_05,
                    data = BIP_merged)

RPS17_residuals <- as.data.frame(RPS17_resid_mod$residuals)

colnames(RPS17_residuals)[1] <- "RPS17_residuals"


BIP_merged$RPS17_residuals <- as.data.frame(scale(RPS17_residuals$RPS17_residuals))

## Residualise the FES and GRIN2A PES

ggplot2.density(data=BIP_merged, 
                xName='RPS17_residuals', groupName='Phenotype', 
                alpha=0.5, fillGroupDensity=TRUE, addMeanLine=TRUE, meanLineColor="black", xtitle="Distribution of residualised RPS17 network PES", 
                ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                densityAlpha=0.3, densityLineSize=4, groupColors=c('darkorchid1','#999999'))



Decile_resid_df <- BIP_merged %>%
  mutate_at(vars(ends_with("residuals")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(ends_with("5")), .funs = list(~ntile(.,10)))

table(Decile_resid_df$RPS17_residuals, Decile_resid_df$Phenotype)

table(Decile_resid_df$RPS17_0_05, Decile_resid_df$Phenotype)

RPS17_decile_discrepancy <- Decile_resid_df %>% filter(RPS17_0_05 == 10 & RPS17_residuals < 10) %>%
  select(Phenotype, RPS17_0_05, RPS17_residuals) %>% filter(Phenotype == "BIP")


## Plot raw FES PES vs residualised FES PES

quantile(BIP_merged$RPS17_0_05, probs = seq(0, 1, 0.1))

quantile(BIP_merged$RPS17_residuals, probs = seq(0, 1, 0.1))


BIP_merged$Plot_var_col <- as.factor(ifelse(BIP_merged$RPS17_0_05 > 1.238206053 & BIP_merged$RPS17_residuals > 1.22916778 , 1, 0))

ggplot(data = BIP_merged,
       aes(x = RPS17_0_05, y=RPS17_residuals, colour = Plot_var_col)) +
  geom_point(alpha = 0.9, size = 1.4) +
  geom_vline(xintercept = 1.238206053 , linetype="longdash") +
  geom_hline(yintercept=1.22916778, linetype="longdash") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Raw RPS17 network PES") +
  ylab("Residualised RPS17 network PES") +
  scale_color_manual(name="Plot_var_col", values=c('#999999','indianred2'))
