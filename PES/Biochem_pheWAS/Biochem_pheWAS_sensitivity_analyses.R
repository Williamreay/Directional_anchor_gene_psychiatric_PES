################################

## Effect of PES/PRS on blood/urine biochemical traits amongst UKBB participants who did not self-report mental illness in the MHQ

## Senstivity analyses

## William Reay (2021)

#################################

library(data.table)
library(dplyr)
library(easyGgplot2)
library(RNOmni)
library(RColorBrewer)
library(ggpubr)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem/")

## Read in df

Merged_biochem_data <- fread("Merged_all_scores_all_biochem.txt")

Merged_biochem_data$Batch <- as.factor(Merged_biochem_data$Batch)

## Read in BMI data

BMI <- fread("BMI_UKBB.txt")

BMI <- rename(BMI, "IID"="f.eid", "BMI"="f.21001.0.0")

## Read in statin data 

Statins <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/Statins_UKBB.tab.txt", header = T)

Statins <- rename(Statins, "IID"="f.eid")

## Merge

Senstivity_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(Merged_biochem_data, BMI, Statins))



## Test adjustment for PRS

BIP_testing <- fread("Biochem_sensitivity_analyses/BIP_list_to_test.txt", header = T)

BIP_testing <- as.list(BIP_testing$Input)

BIP_scores <- fread("Biochem_sensitivity_analyses/BIP_scores_to_test.txt", header = T)

BIP_scores <- as.list(BIP_scores$Input)

## BIP

BIP_PRS_adjustment <- function(biochem_col, score_col, df) {
  biochem_df <- df[!is.na(df[[biochem_col]]),]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  fmla <- as.formula(paste(biochem_col, "~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + scaled_score + BIP_PRS +  Batch"))
  mod <- lm(fmla, data = biochem_df)
  return(summary(mod))
}

Run_BIP <- mapply(BIP_PRS_adjustment, biochem_col = BIP_testing, score_col=BIP_scores, 
                  MoreArgs = list(Senstivity_merged),
                  SIMPLIFY = TRUE)

Run_BIP_extract <- apply(Run_BIP, 2, function(x) return(as.data.frame(x$coefficients)[15, 1:4]))

## SZ

SZ_testing <- fread("Biochem_sensitivity_analyses/SZ_list_to_test.txt", header = T)

SZ_testing <- as.list(SZ_testing$Input)

SZ_PRS_adjustment <- function(biochem_col, score_col, df) {
  biochem_df <- df[!is.na(df[[biochem_col]]),]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  fmla <- as.formula(paste(biochem_col, "~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + scaled_score + SZ_PRS +  Batch"))
  mod <- lm(fmla, data = biochem_df)
  return(summary(mod))
}

Run_SZ <- sapply(SZ_testing, SZ_PRS_adjustment,  score_col="SZ_FADS1_PES", df=Senstivity_merged)

Run_SZ_extract <- apply(Run_SZ, 2, function(x) return(as.data.frame(x$coefficients)[15, 1:4]))

## Test natural log transformation of outcome variable

All_scores <- fread("Biochem_sensitivity_analyses/All_sensitivity_scores.txt", header = T)

All_scores <- as.list(All_scores$Input)

All_testing <- fread("Biochem_sensitivity_analyses/All_sensitivity_biochem.txt", header = T)

All_testing <- as.list(All_testing$Input)

Ln_transform <- function(biochem_col, score_col, df) {
  biochem_df <- df[!is.na(df[[biochem_col]]),]
  biochem_df$ln_biochem <- log(biochem_df[[biochem_col]])
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  fmla <- as.formula("ln_biochem ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + scaled_score +  Batch")
  mod <- lm(fmla, data = biochem_df)
  return(summary(mod))
}

Run_ln <- mapply(Ln_transform, biochem_col = All_testing, score_col=All_scores, 
                  MoreArgs = list(Senstivity_merged),
                  SIMPLIFY = TRUE)

Run_ln_extract <- apply(Run_ln, 2, function(x) return(as.data.frame(x$coefficients)[15, 1:4]))

Output <- data.frame()
for (i in 1:length(Run_ln_extract)) {
  Output <- rbind(Output, Run_ln_extract[[i]])
}

write.table(Output, file="Biochem_sensitivity_analyses/Ln_transform_results.txt",
            sep = "\t", row.names = F, quote = F)


## Test effect for statin usage

Statin_adjust <- function(biochem_col, score_col, df) {
  biochem_df <- df[!is.na(df[[biochem_col]]),]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  fmla <- as.formula(paste(biochem_col," ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + scaled_score + Statins + Batch"))
  mod <- lm(fmla, data = biochem_df)
  return(summary(mod))
}

Run_statin <- mapply(Statin_adjust, biochem_col = All_testing, score_col=All_scores, 
                 MoreArgs = list(Senstivity_merged),
                 SIMPLIFY = TRUE)

Run_statin_extract <- apply(Run_statin, 2, function(x) return(as.data.frame(x$coefficients)[15, 1:4]))

S_output <- data.frame()
for (i in 1:length(Run_statin_extract)) {
  S_output <- rbind(S_output, Run_statin_extract[[i]])
}

write.table(S_output, file="Biochem_sensitivity_analyses/Statin_covariation_results.txt",
            sep = "\t", row.names = F, quote = F)

## Inverse-rank normal transformation

IRNT <- function(biochem_col, score_col, df) {
  biochem_df <- df[!is.na(df[[biochem_col]]),]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  fmla <- as.formula(paste(biochem_col," ~ Sex*Age + Age2"))
  mod <- lm(fmla, data = biochem_df)
  Biochem_resid <- as.numeric(mod$residuals)
  ## Blom transformation
  Biochem_residuals_NORMALISED <- as.data.frame(RankNorm(Biochem_resid), k=3/8)
  biochem_df$Norm_resid <- Biochem_residuals_NORMALISED$'RankNorm(Biochem_resid)'
  Resid_mod <- lm(Norm_resid ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + scaled_score + Batch, data = biochem_df)
  return(summary(Resid_mod))
}

Run_IRNT <- mapply(IRNT, biochem_col = All_testing, score_col=All_scores, 
                     MoreArgs = list(Senstivity_merged),
                     SIMPLIFY = TRUE)

Run_IRNT_extract <- apply(Run_IRNT, 2, function(x) return(as.data.frame(x$coefficients)[12, 1:4]))

IRNT_output <- data.frame()
for (i in 1:length(Run_IRNT_extract)) {
  IRNT_output <- rbind(IRNT_output, Run_IRNT_extract[[i]])
}

IRNT_output$Biochem <- unlist(All_testing)
IRNT_output$Score <- unlist(All_scores)

write.table(IRNT_output, file="Biochem_sensitivity_analyses/IRNT_results.txt",
            sep = "\t", row.names = F, quote = F)

## Make barplot of deciles for HDL with each of the scores

HDL_df <- Senstivity_merged %>% filter(!is.na(f.30760.0.0))

HDL_df <- HDL_df %>% select(f.30760.0.0, SZ_FADS1_PES, BIP_FADS1_PES,
                            SZ_PRS, BIP_PRS, Sex, Age2, Age)

## Scale scores

HDL_df[,c(2:5)] <- lapply(HDL_df[,c(2:5)], function(x) c(scale(x)))

## Make deciles of scores 

HDL_df <- HDL_df %>%
  mutate_at(vars(starts_with("SZ")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("BIP")), .funs = list(~ntile(.,10)))


HDL_df <- rename(HDL_df, "HDL"="HDL_IRNT")

HDL_df <- HDL_df %>% filter(SZ_FADS1_PES == 1 | SZ_FADS1_PES == 5 |
                              SZ_FADS1_PES == 10)

HDL_df <- HDL_df %>% filter(BIP_FADS1_PES == 1 | BIP_FADS1_PES == 5 |
                              BIP_FADS1_PES == 10)

HDL_df <- HDL_df %>% filter(SZ_PRS == 1 | SZ_PRS == 5 |
                              SZ_PRS == 10)

HDL_df <- HDL_df %>% filter(BIP_PRS == 1 | BIP_PRS == 5 |
                              BIP_PRS == 10)

HDL_df$SZ_FADS1_PES <- as.factor(HDL_df$SZ_FADS1_PES)
HDL_df$BIP_FADS1_PES <- as.factor(HDL_df$BIP_FADS1_PES)
HDL_df$SZ_PRS <- as.factor(HDL_df$SZ_PRS)
HDL_df$BIP_PRS <- as.factor(HDL_df$BIP_PRS)


## Make plots

BIP_HDL <- ggplot(HDL_df, aes(x=BIP_FADS1_PES, y=f.30760.0.0, fill=BIP_FADS1_PES)) +
  geom_boxplot() +
  theme_bw() +
  xlab("") +
  ylab("HDL (mmol/L))") +
  labs(fill = "Score name") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1.478, linetype = "longdash") +
  scale_x_discrete(labels=c("1" = "≤ 10th percentile", "5" = "50th-60th percentile",
                            "10" = "≥ 90th percentile")) +
  scale_fill_brewer(name = "Blues") +
  ggtitle("BIP FADS1 PES - HDL")

SZ_HDL <- ggplot(HDL_df, aes(x=SZ_FADS1_PES, y=f.30760.0.0,
                             fill = SZ_FADS1_PES)) +
  geom_boxplot() +
  theme_bw() +
  xlab("") +
  ylab("HDL (mmol/L)") +
  labs(fill = "Score name") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1.478, linetype = "longdash") +
  scale_x_discrete(labels=c("1" = "≤ 10th percentile", "5" = "50th-60th percentile",
                            "10" = "≥ 90th percentile")) +
  scale_fill_brewer(palette = "Blues") +
  ggtitle("SZ FADS1 PES - HDL")


HDL_all <- ggarrange(BIP_HDL, SZ_HDL)


