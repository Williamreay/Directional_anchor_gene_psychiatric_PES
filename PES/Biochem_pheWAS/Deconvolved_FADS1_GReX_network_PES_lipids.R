################################

## Effect of deconvolved FADS1 PES and FADS1 GREx on lipids

## William Reay (2021)

#################################

library(data.table)
library(dplyr)
library(easyGgplot2)
library(RColorBrewer)
library(ggpubr)
library(readxl)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/")

## Define concatenation function

PES_concat <- function(df, name){
  ## Format df with unique colnames
  Colnames_df <- make.unique(names(df))
  colnames(df) <- Colnames_df
  ## Select columns for PES or PRS
  Score_df <- df %>% select(IID, contains("AVG"))
  ## Sum chromosome-wise score
  Score_df[[name]] <- rowSums(Score_df[,-1])
  Score_df <- Score_df %>% select(IID, name)
  return(Score_df)
}

## SZ network PES (no FADS1) ##

SZ_FADS1 <- fread("BIP_best_con_or_lib/FADS1_GRex_backbone_deconvolved/Scores/All_chr_SZ_FADS1_network_no_FADS1_region.tab.txt", header = T)

SZ_FADS1_concat <- PES_concat(SZ_FADS1, "SZ_FADS1")

write.table(SZ_FADS1_concat, file="BIP_best_con_or_lib/FADS1_GRex_backbone_deconvolved/Scores/Concatenated_SZ_FADS1_no_FADS1.txt",
            row.names = F, quote = F, sep = "\t")

## BIP network PES (no FADS1) ##

BIP_FADS1 <- fread("BIP_best_con_or_lib/FADS1_GRex_backbone_deconvolved/Scores/All_chr_BIP_FADS1_network_no_FADS1_region.tab.txt", header = T)

BIP_FADS1_concat <- PES_concat(BIP_FADS1, "BIP_FADS1")

write.table(BIP_FADS1_concat, file="BIP_best_con_or_lib/FADS1_GRex_backbone_deconvolved/Scores/Concatenated_BIP_FADS1_no_FADS1.txt",
            row.names = F, quote = F, sep = "\t")

## FADS1 GReX ##

FADS1_GReX <- fread("BIP_best_con_or_lib/FADS1_GRex_backbone_deconvolved/Scores/UKBB_FADS1_GReX_11.sscore", header = T)

## Read in df

Merged_biochem_data <- fread("UKBB_MHQ_no_self_reported_biochem/Merged_all_scores_all_biochem.txt")

Merged_biochem_data$Batch <- as.factor(Merged_biochem_data$Batch)

## Merge all 

FADS1_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                            list(Merged_biochem_data, SZ_FADS1_concat, BIP_FADS1_concat, FADS1_GReX))

## HDL ##

## Make df with non-missing HDL

HDL_df <- FADS1_merged %>% filter(!is.na(f.30760.0.0))

HDL_df <- HDL_df %>%
  mutate_at(vars(ends_with("FADS1")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("FADS1")), .funs = list(~ntile(.,10)))

HDL_df$Low_GReX <- ifelse(HDL_df$FADS1_GReX_AVG == 1, 1, 0)
HDL_df$High_PES <- ifelse(HDL_df$SZ_FADS1 == 10, 1, 0)
HDL_df$High_BIP_PES <- ifelse(HDL_df$BIP_FADS1 == 10, 1, 0)
HDL_df$High_PES_low_GReX <- ifelse(HDL_df$FADS1_GReX_AVG == 1 & HDL_df$SZ_FADS1 == 1, 1, 0)

HDL_df$High_BIP_PES_low_GReX <- ifelse(HDL_df$FADS1_GReX_AVG == 1 & HDL_df$BIP_FADS1 == 1, 1, 0)

HDL_SZ_GReX <- lm(f.30760.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                 PC8 + PC9 + PC10 + High_PES_low_GReX +  Batch, data = HDL_df)

HDL_GReX <- lm(f.30760.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + Low_GReX +  Batch, data = HDL_df)

HDL_SZ <- lm(f.30760.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + High_PES +  Batch, data = HDL_df)

HDL_BIP_GReX_PES <- lm(f.30760.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                 PC8 + PC9 + PC10 + High_BIP_PES_low_GReX +  Batch, data = HDL_df)

## Add full PES

HDL_GReX_mod <- lm(f.30760.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                 PC8 + PC9 + PC10 + FADS1_GReX_AVG + BIP_FADS1_PES + Batch, data = HDL_df)

anova(HDL_GReX_mod, test="F")

## Make FP

FP_input <- read_excel("BIP_best_con_or_lib/FADS1_GRex_backbone_deconvolved/GReX_FP.xlsx")

FP_input$U_CI <- FP_input$Beta + (1.96*FP_input$SE)

FP_input$L_CI <- FP_input$Beta - (1.96*FP_input$SE)

ggplot(data = FP_input, aes(x=Score, y=Beta, ymin=U_CI, ymax=L_CI, colour=Score)) + geom_pointrange() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  geom_hline(yintercept=0, lty=2) +
  ylab("Effect on HDL (mmol/L) per SD in score") +
  scale_color_brewer(palette="Paired")

HDL_PES_df <- HDL_df %>% filter(SZ_FADS1 == 1 | SZ_FADS1 == 5 |
                                SZ_FADS1 == 10)

HDL_FADS1_df <- HDL_df %>% filter(FADS1_GReX_AVG == 1 |
                                  FADS1_GReX_AVG == 5 | FADS1_GReX_AVG == 10)



HDL_PES_df$SZ_FADS1 <- as.factor(HDL_PES_df$SZ_FADS1)
HDL_FADS1_df$FADS1_GReX_AVG <- as.factor(HDL_FADS1_df$FADS1_GReX_AVG)


## Scale scores

HDL_df$SZ_FADS1 <- as.numeric(scale(HDL_df$SZ_FADS1))
HDL_df$BIP_FADS1 <- as.numeric(scale(HDL_df$BIP_FADS1))
HDL_df$FADS1_GReX_AVG <- as.numeric(scale(HDL_df$FADS1_GReX_AVG))

## Test association between GReX and HDL

HDL_GReX <- lm(f.30760.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
  PC8 + PC9 + PC10 + FADS1_GReX_AVG +  Batch, data = HDL_df)

## Test association between networks

HDL_SZ_network <- lm(f.30760.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                 PC8 + PC9 + PC10 + SZ_FADS1 + Batch, data = HDL_df)

HDL_BIP_network <- lm(f.30760.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                       PC8 + PC9 + PC10 + BIP_FADS1 +  Batch, data = HDL_df)

## Test for interaction effects

GReX_int_HDL <- lm(Norm_resid ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                     PC8 + PC9 + PC10 + FADS1_GReX_AVG*SZ_FADS1 +  Batch, data = HDL_df)

## Triglycerides ##

## Make df with non-missing HDL

TG_df <- FADS1_merged %>% filter(!is.na(f.30870.0.0))

## Identify individuals with elevated PES but low GReX

TG_df <- TG_df %>%
  mutate_at(vars(ends_with("FADS1")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("FADS1")), .funs = list(~ntile(.,10)))

TG_df$Low_GReX <- ifelse(TG_df$FADS1_GReX_AVG == 1, 1, 0)
TG_df$High_PES <- ifelse(TG_df$BIP_FADS1 == 10, 1, 0)

TG_BIP_decile <- lm(f.30870.0.0 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                      PC8 + PC9 + PC10 + High_PES +  Batch, data = TG_df)

TG_GReX_decile <- lm(f.30870.0.0 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                       PC8 + PC9 + PC10 + Low_GReX +  Batch, data = TG_df)


## Individuals with high PES and low GReX

TG_df$High_PES_and_low_GReX <- ifelse(TG_df$Low_GReX == 1 & TG_df$High_PES == 1, 1, 0)

TG_subset_df <- TG_df %>% filter(TG_df$High_PES_and_low_GReX == 1 | TG_df$Low_GReX == 1)

TG_both_PES_GReX <- lm(f.30870.0.0  ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                         PC8 + PC9 + PC10 + High_PES_and_low_GReX +  Batch, data = TG_subset_df)

## IRNT to boost power

TG_mod <- lm(f.30870.0.0 ~ Sex*Age + Age2, data = TG_df)

TG_resid <- as.numeric(TG_mod$residuals)

TG_residuals_NORMALISED <- as.data.frame(RankNorm(TG_resid), k=3/8)
TG_df$Norm_resid <- TG_residuals_NORMALISED$'RankNorm(TG_resid)'


## Scale scores

TG_df$SZ_FADS1 <- as.numeric(scale(TG_df$SZ_FADS1))
TG_df$BIP_FADS1 <- as.numeric(scale(TG_df$BIP_FADS1))
TG_df$FADS1_GReX_AVG <- as.numeric(scale(TG_df$FADS1_GReX_AVG))

## Test association between GReX and HDL

TG_GReX <- lm(Norm_resid ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                 PC8 + PC9 + PC10 + FADS1_GReX_AVG +  Batch, data = TG_df)

## Test association between networks

TG_SZ_network <- lm(Norm_resid ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                       PC8 + PC9 + PC10 + SZ_FADS1 +  Batch, data = TG_df)

TG_BIP_network <- lm(Norm_resid ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                        PC8 + PC9 + PC10 + BIP_FADS1 +  Batch, data = TG_df)

## Test for interaction effects

GReX_int_TG <- lm(Norm_resid ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                     PC8 + PC9 + PC10 + FADS1_GReX_AVG*BIP_FADS1 +  Batch, data = TG_df)


TG_PES_df <- TG_df %>% filter(BIP_FADS1 == 1 | BIP_FADS1 == 5 |
                              BIP_FADS1 == 10)

TG_FADS1_df <- TG_df %>% filter(FADS1_GReX_AVG == 1 |
  FADS1_GReX_AVG == 5 | FADS1_GReX_AVG == 10)


TG_PES_df$BIP_FADS1 <- as.factor(TG_PES_df$BIP_FADS1)
TG_FADS1_df$FADS1_GReX_AVG <- as.factor(TG_FADS1_df$FADS1_GReX_AVG)

BIP_TG <- ggplot(TG_PES_df, aes(x=BIP_FADS1, y=Norm_resid, fill=BIP_FADS1)) +
  geom_boxplot() +
  theme_bw() +
  xlab("") +
  ylab("TG (SD)") +
  labs(fill = "Score name") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "longdash") +
  scale_x_discrete(labels=c("1" = "≤ 10th percentile", "5" = "50th-60th percentile",
                            "10" = "≥ 90th percentile")) +
  scale_fill_brewer(name = "Blues") +
  ggtitle("BIP FADS1 PES - TG")

BIP_GReX <- ggplot(TG_FADS1_df, aes(x=FADS1_GReX_AVG, y=f.30870.0.0, fill=FADS1_GReX_AVG)) +
  geom_boxplot() +
  theme_bw() +
  xlab("") +
  ylab("TG (mmol/L)") +
  labs(fill = "Score name") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "longdash") +
  scale_x_discrete(labels=c("1" = "≤ 10th percentile", "5" = "50th-60th percentile",
                            "10" = "≥ 90th percentile")) +
  scale_fill_brewer(name = "Blues") +
  ggtitle("BIP FADS1 GReX - TG")


## Total cholesterol ##

## Make df with non-missing HDL

TC_df <- FADS1_merged %>% filter(!is.na(f.30690.0.0))

## Scale scores

TC_df$SZ_FADS1 <- as.numeric(scale(TC_df$SZ_FADS1))
TC_df$BIP_FADS1 <- as.numeric(scale(TC_df$BIP_FADS1))
TC_df$FADS1_GReX_AVG <- as.numeric(scale(TC_df$FADS1_GReX_AVG))

## Test association between GReX and HDL

TC_GReX <- lm(f.30690.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                PC8 + PC9 + PC10 + FADS1_GReX_AVG +  Batch, data = HDL_df)

## Test association between networks

TC_SZ_network <- lm(f.30690.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                      PC8 + PC9 + PC10 + SZ_FADS1 +  Batch, data = HDL_df)

TC_BIP_network <- lm(f.30690.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                       PC8 + PC9 + PC10 + BIP_FADS1 +  Batch, data = HDL_df)

## Test for interaction effects

GReX_int_TC <- lm(f.30690.0.0 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + FADS1_GReX_AVG*BIP_FADS1 +  Batch, data = HDL_df)
