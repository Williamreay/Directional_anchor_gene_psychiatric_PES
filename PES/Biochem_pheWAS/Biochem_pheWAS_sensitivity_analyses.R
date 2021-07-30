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

## Read in statin data 

Statins <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/Statins_UKBB.tab.txt", header = T)

Statins <- rename(Statins, "IID"="f.eid")

## Merge

Senstivity_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(Merged_biochem_data, Statins))


## SZ ##

## Test adjustment for PRS

SZ_testing <- fread("Biochem_sensitivity_analyses/SZ_list_to_test.txt", header = T)

SZ_testing <- as.list(SZ_testing$Input)

SZ_scores <- fread("Biochem_sensitivity_analyses/SZ_scores_to_test.txt", header = T)

SZ_scores <- as.list(SZ_scores$Input)

SZ_PRS_adjustment <- function(biochem_col, score_col, df) {
  biochem_df <- df[!is.na(df[[biochem_col]]),]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  fmla <- as.formula(paste(biochem_col, "~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + scaled_score + SZ_PRS +  Batch"))
  mod <- lm(fmla, data = biochem_df)
  return(summary(mod))
}

Run_SZ <- mapply(SZ_PRS_adjustment, biochem_col = SZ_testing, score_col=SZ_scores, 
                  MoreArgs = list(Senstivity_merged),
                  SIMPLIFY = TRUE)

Run_SZ_extract <- apply(Run_SZ, 2, function(x) return(as.data.frame(x$coefficients)[15, 1:4]))

Run_SZ_results <- data.frame()

for (i in 1:length(Run_SZ_extract)) {
  Run_SZ_results <- rbind(Run_SZ_results, Run_SZ_extract[[i]])
}
Run_SZ_results$Pheno <- unlist(SZ_testing)
Run_SZ_results$Score <- unlist(SZ_scores)

write.csv(Run_SZ_results, file="Biochem_sensitivity_analyses/PRS_adjustment_SZ_scores_FDR_sig.csv",
          row.names = F, quote = F)

## Interaction model

SZ_PRS_interaction <- function(biochem_col, score_col, df) {
  biochem_df <- df[!is.na(df[[biochem_col]]),]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  fmla <- as.formula(paste(biochem_col, "~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + scaled_score*SZ_PRS +  Batch"))
  mod <- lm(fmla, data = biochem_df)
  return(summary(mod))
}

Run_SZ_int <- mapply(SZ_PRS_interaction, biochem_col = SZ_testing, score_col=SZ_scores, 
                  MoreArgs = list(Senstivity_merged),
                  SIMPLIFY = TRUE)

Run_SZ_extract_int <- apply(Run_SZ_int, 2, function(x) return(as.data.frame(x$coefficients)[123, 1:4]))

## BIP

BIP_testing <- fread("Biochem_sensitivity_analyses/BIP_list_to_test.txt", header = T)

BIP_testing <- as.list(BIP_testing$Input)

BIP_PRS_adjustment <- function(biochem_col, score_col, df) {
  biochem_df <- df[!is.na(df[[biochem_col]]),]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  fmla <- as.formula(paste(biochem_col, "~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + scaled_score + BIP_PRS +  Batch"))
  mod <- lm(fmla, data = biochem_df)
  return(summary(mod))
}

Run_BIP <- sapply(BIP_testing, BIP_PRS_adjustment,  score_col="BIP_FADS1_PES", df=Senstivity_merged)

Run_BIP_extract <- apply(Run_BIP, 2, function(x) return(as.data.frame(x$coefficients)[15, 1:4]))

Run_BIP_results <- data.frame()

for (i in 1:length(Run_BIP_extract)) {
  Run_BIP_results <- rbind(Run_BIP_results, Run_BIP_extract[[i]])
}
Run_BIP_results$Pheno <- unlist(BIP_testing)
Run_BIP_results$Score <- unlist(BIP_scores)

write.csv(Run_BIP_results, file="Biochem_sensitivity_analyses/PRS_adjustment_BIP_scores_FDR_sig.csv",
          row.names = F, quote = F)

## Interaction term

BIP_PRS_int <- function(biochem_col, score_col, df) {
  biochem_df <- df[!is.na(df[[biochem_col]]),]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  fmla <- as.formula(paste(biochem_col, "~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + BIP_PRS*scaled_score +  Batch"))
  mod <- lm(fmla, data = biochem_df)
  return(summary(mod))
}

Run_BIP_int <- sapply(BIP_testing, BIP_PRS_int,  score_col="BIP_FADS1_PES", df=Senstivity_merged)

Run_BIP_extract_int <- apply(Run_BIP_int, 2, function(x) return(as.data.frame(x$coefficients)[123, 1:4]))

## Test natural log transformation of outcome variable

All_scores <- fread("Biochem_sensitivity_analyses/All_sensitvity_scores.txt", header = T)

All_scores <- as.list(All_scores$Input)

All_testing <- fread("Biochem_sensitivity_analyses/All_sensitvity_biochem.txt", header = T)

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

Output$Pheno <- unlist(All_testing)
Output$Score <- unlist(All_scores)

write.csv(Output, file="Biochem_sensitivity_analyses/Ln_transform_results.csv",
            row.names = F, quote = F)


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

S_output$Pheno <- unlist(All_testing)
S_output$Score <- unlist(All_scores)

write.csv(S_output, file="Biochem_sensitivity_analyses/Statin_covariation_results.csv",
            row.names = F, quote = F)

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

write.csv(IRNT_output, file="Biochem_sensitivity_analyses/IRNT_results.csv",
            row.names = F, quote = F)

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


HDL_df_SZ_PES <- HDL_df %>% filter(SZ_FADS1_PES == 1 | SZ_FADS1_PES == 5 |
                              SZ_FADS1_PES == 10)

HDL_df_BIP_PES <- HDL_df %>% filter(BIP_FADS1_PES == 1 | BIP_FADS1_PES == 5 |
                              BIP_FADS1_PES == 10)

HDL_df_SZ_PRS <- HDL_df %>% filter(SZ_PRS == 1 | SZ_PRS == 5 |
                              SZ_PRS == 10)

HDL_df_BIP_PRS <- HDL_df %>% filter(BIP_PRS == 1 | BIP_PRS == 5 |
                              BIP_PRS == 10)

HDL_df_SZ_PES$SZ_FADS1_PES <- as.factor(HDL_df_SZ_PES$SZ_FADS1_PES)
HDL_df_BIP_PES$BIP_FADS1_PES <- as.factor(HDL_df_BIP_PES$BIP_FADS1_PES)
HDL_df_SZ_PRS$SZ_PRS <- as.factor(HDL_df_SZ_PRS$SZ_PRS)
HDL_df_BIP_PRS$BIP_PRS <- as.factor(HDL_df_BIP_PRS$BIP_PRS)


## Make plots

BIP_HDL <- ggplot(HDL_df_BIP_PES, aes(x=BIP_FADS1_PES, y=f.30760.0.0, fill=BIP_FADS1_PES)) +
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

SZ_HDL <- ggplot(HDL_df_SZ_PES, aes(x=SZ_FADS1_PES, y=f.30760.0.0,
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

## FES i IGF1 ##

IGF1_df <- Senstivity_merged %>% filter(!is.na(f.30770.0.0))

IGF1_df <- IGF1_df %>% select(f.30770.0.0, SZ_FES_PES, BIP_FES_PES,
                              SZ_PRS, BIP_PRS, Sex, Age2, Age)

## Scale scores

IGF1_df[,c(2:5)] <- lapply(IGF1_df[,c(2:5)], function(x) c(scale(x)))

## Make deciles of scores 

IGF1_df <- IGF1_df %>%
  mutate_at(vars(starts_with("SZ")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("BIP")), .funs = list(~ntile(.,10)))


IGF1_df_SZ_PES <- IGF1_df %>% filter(SZ_FES_PES == 1 | SZ_FES_PES == 5 |
                                       SZ_FES_PES == 10)

IGF1_df_BIP_PES <- IGF1_df %>% filter(BIP_FES_PES == 1 | BIP_FES_PES == 5 |
                                        BIP_FES_PES == 10)

IGF1_df_SZ_PRS <- IGF1_df %>% filter(SZ_PRS == 1 | SZ_PRS == 5 |
                                       SZ_PRS == 10)

IGF1_df_BIP_PRS <- IGF1_df %>% filter(BIP_PRS == 1 | BIP_PRS == 5 |
                                        BIP_PRS == 10)

IGF1_df_SZ_PES$SZ_FES_PES <- as.factor(IGF1_df_SZ_PES$SZ_FES_PES)
IGF1_df_BIP_PES$BIP_FES_PES <- as.factor(IGF1_df_BIP_PES$BIP_FES_PES)
IGF1_df_SZ_PRS$SZ_PRS <- as.factor(IGF1_df_SZ_PRS$SZ_PRS)
IGF1_df_BIP_PRS$BIP_PRS <- as.factor(IGF1_df_BIP_PRS$BIP_PRS)


## Make plots

BIP_IGF1 <- ggplot(IGF1_df_BIP_PES, aes(x=BIP_FES_PES, y=f.30770.0.0, fill=BIP_FES_PES)) +
  geom_boxplot() +
  theme_bw() +
  xlab("") +
  ylab("IGF1 (nmol/L))") +
  labs(fill = "Score name") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 21.88, linetype = "longdash") +
  scale_x_discrete(labels=c("1" = "≤ 10th percentile", "5" = "50th-60th percentile",
                            "10" = "≥ 90th percentile")) +
  scale_fill_brewer(name = "Blues") +
  ggtitle("BIP FES PES - IGF1")

SZ_IGF1 <- ggplot(IGF1_df_SZ_PES, aes(x=SZ_FES_PES, y=f.30770.0.0,
                                      fill = SZ_FES_PES)) +
  geom_boxplot() +
  theme_bw() +
  xlab("") +
  ylab("IGF1 (nmol/L)") +
  labs(fill = "Score name") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 21.88, linetype = "longdash") +
  scale_x_discrete(labels=c("1" = "≤ 10th percentile", "5" = "50th-60th percentile",
                            "10" = "≥ 90th percentile")) +
  scale_fill_brewer(palette = "Blues") +
  ggtitle("SZ FES PES - IGF1")


IGF1_all <- ggarrange(BIP_IGF1, SZ_IGF1)

