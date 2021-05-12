#############################################

## Association of PES and PRS with SZ in the ASRB cohort (sumstats with ASRB removed)

## William Reay (2021)

#############################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/ASRB_SZ/")

library(data.table)
library(readxl)
library(dplyr)
library(rcompanion)
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)


## Read in concatenated PES/PRS file

All_SZ_scores <- fread("Combined_ASRB_SZ_PES_PRS.txt", header = T, sep = " ")

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores

## Read in phenotype and covariate data

SZ_pheno <- read_excel("ASRB_case_vs_control.xlsx")

## Merge case vs control df

ASRB_case_control <- merge(All_SZ_scores, SZ_pheno, by = "IID")

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

ASRB_case_control[,c(3, 7, 11, 15, 19, 23:26, 29)]

ASRB_case_control[,c(3, 7, 11, 15, 19, 23:26, 29)] <- lapply(ASRB_case_control[,c(3, 7, 11, 15, 19, 23:26, 29)], function(x) c(scale(x)))

## Merge with CD vs CS df

GoM <- read.csv("CD_vs_CS_GoM1.csv", header = T)

ASRB_GoM <- merge(All_SZ_scores, GoM, by = "IID")

ASRB_GoM[,c(3, 7, 11, 15, 19, 23:26, 29)] <- lapply(ASRB_GoM[,c(3, 7, 11, 15, 19, 23:26, 29)], function(x) c(scale(x)))

## Merge with clozapine df

Cloz <- read.table("Clozapine_ASRB_new_df.txt", header = T)

Cloz_merge <- merge(Cloz, SZ_pheno, by = "IID")

ASRB_cloz <- merge(All_SZ_scores, Cloz_merge, by = "IID")

ASRB_cloz[,c(3, 7, 11, 15, 19, 23:26, 29)] <- lapply(ASRB_cloz[,c(3, 7, 11, 15, 19, 23:26, 29)], function(x) c(scale(x)))

## Define function to test each score by itself and extract beta, SE, and P

## List of sets to test

Scores_to_test <- as.list(colnames(ASRB_case_control[,c(3, 7, 11, 15, 19, 23:26, 29)]))

## Function

PES_PRS_SZ_ASRB <- function(v){
  Score_model <- glm(glue::glue('PHENOTYPE ~ SEX + PC1 + PC2 + PC3 + {v}'), family = "binomial", data = ASRB_case_control)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

SZ_score <- sapply(Scores_to_test, PES_PRS_SZ_ASRB)

## Extract beta, se, z, and p value for each gene

SZ_extract <- apply(SZ_score, 2, function(x) return(as.data.frame(x$coefficients)[6,1:4]))

SZ_results <- data.frame()

for (i in 1:length(SZ_extract)) {
  SZ_results <- rbind(SZ_results, SZ_extract[[i]])
}
rownames(SZ_results) <- Scores_to_test

write.csv(SZ_results, file="Results/SZ_ASRB_PES_PRS_results.csv", quote = F, row.names = T)


## Make forest plot of effect size per SD

SZ_results$Scores <- Scores_to_test


Plot_df <-  SZ_results %>% filter(Scores == "FES_0_5" |
                                    Scores == "RPS17_1" | Scores == "GRIN2A_0_005" |
                                    Scores == "CACNA1C_0_005" | Scores == "PCCB_0_05" |
                                    Scores == "FADS1_0_05")

Plot_df$Score_name <- c("CACNA1C netowrk PES", "FADS1 network PES", 
                        "FES network PES", "GRIN2A network PES",
                        "PCCB network PES", "RPS17 network PES")

## Plot results

Plot_df$OR <- exp(Plot_df$Estimate)
Plot_df$UOR <- Plot_df$OR + 1.96*Plot_df$`Std. Error`
Plot_df$LOR <- Plot_df$OR - 1.96*Plot_df$`Std. Error`

FP_SZ <- ggplot(data = Plot_df, aes(x=Score_name, y=OR, ymin=LOR, ymax=UOR, colour=Score_name)) +
  geom_pointrange() +
  geom_hline(yintercept=1, lty=2) +
  coord_flip() +
  ylab("Schizophrenia OR [95% CI] per SD increase in PES") +
  theme_bw() +
  theme(legend.position = "null", axis.title.y = element_blank())

FP_SZ + scale_color_brewer(palette="Dark2")


## GoM function

## Function

PES_PRS_SZ_GoM <- function(v){
  Score_model <- glm(glue::glue('GoM_coded ~ SEX + PC1 + PC2 + PC3 + {v}'), family = "binomial", data = ASRB_GoM)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

GoM_score <- sapply(Scores_to_test, PES_PRS_SZ_GoM)

## Extract beta, se, z, and p value for each gene

GoM_extract <- apply(GoM_score, 2, function(x) return(as.data.frame(x$coefficients)[6,1:4]))

GoM_results <- data.frame()

for (i in 1:length(GoM_extract)) {
  GoM_results <- rbind(GoM_results, GoM_extract[[i]])
}
rownames(GoM_results) <- Scores_to_test

write.csv(GoM_results, file="Results/SZ_ASRB_CD_CS_PES_PRS_results.csv", quote = F, row.names = T)

## Check covariation of RPS17 score

RPS17_PES_CD_null <- glm(GoM_coded ~ SEX + PC1 + PC2 + PC3 + PRS_1,
                         family = "binomial", data = ASRB_GoM)

RPS17_PES_CD_full <- glm(GoM_coded ~ SEX + PC1 + PC2 + PC3 + PRS_1 + RPS17_1,
                         family = "binomial", data = ASRB_GoM)

anova(RPS17_PES_CD_full, RPS17_PES_CD_null, test = "Chisq")


## Make forest plot of effect size per SD

GoM_results$Scores <- Scores_to_test


GoM_Plot_df <-  GoM_results %>% filter(Scores == "FES_0_5" |
                                    Scores == "RPS17_1" | Scores == "GRIN2A_0_005" |
                                    Scores == "CACNA1C_0_005" | Scores == "PCCB_0_05" |
                                    Scores == "FADS1_0_05")

GoM_Plot_df$Score_name <- c("CACNA1C netowrk PES", "FADS1 network PES", 
                        "FES network PES", "GRIN2A network PES",
                        "PCCB network PES", "RPS17 network PES")

## Plot results

GoM_Plot_df$OR <- exp(GoM_Plot_df$Estimate)
GoM_Plot_df$UOR <- GoM_Plot_df$OR + 1.96*GoM_Plot_df$`Std. Error`
GoM_Plot_df$LOR <- GoM_Plot_df$OR - 1.96*GoM_Plot_df$`Std. Error`

FP_GoM <- ggplot(data = GoM_Plot_df, aes(x=Score_name, y=OR, ymin=LOR, ymax=UOR, colour=Score_name)) +
  geom_pointrange() +
  geom_hline(yintercept=1, lty=2) +
  coord_flip() +
  ylab("Cognitive deficit OR [95% CI] per SD increase in PES") +
  theme_bw() +
  theme(legend.position = "null", axis.title.y = element_blank())

FP_GoM + scale_color_brewer(palette="Paired")


## Clozapine function

## Function

PES_PRS_SZ_Cloz <- function(v){
  Score_model <- glm(glue::glue('Clozapine ~ SEX + PC1 + PC2 + PC3 + {v}'), family = "binomial", data = ASRB_cloz)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

Cloz_score <- sapply(Scores_to_test, PES_PRS_SZ_Cloz)

## Extract beta, se, z, and p value for each gene

Cloz_extract <- apply(Cloz_score, 2, function(x) return(as.data.frame(x$coefficients)[6,1:4]))

Cloz_results <- data.frame()

for (i in 1:length(Cloz_extract)) {
  Cloz_results <- rbind(Cloz_results, Cloz_extract[[i]])
}
rownames(Cloz_results) <- Scores_to_test

write.csv(Cloz_results, file="Results/SZ_ASRB_clozapine_PES_PRS_results.csv", quote = F, row.names = T)


## Onset age 

## Function

PES_PRS_SZ_onset <- function(v){
  Score_model <- glm(glue::glue('Coded_onset ~ SEX + PC1 + PC2 + PC3 + {v}'), family = "binomial", data = ASRB_cloz_onset)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

Onset_score <- sapply(Scores_to_test, PES_PRS_SZ_onset)

## Extract beta, se, z, and p value for each gene

Onset_extract <- apply(Onset_score, 2, function(x) return(as.data.frame(x$coefficients)[6,1:4]))

Onset_results <- data.frame()

for (i in 1:length(Onset_extract)) {
  Onset_results <- rbind(Onset_results, Onset_extract[[i]])
}
rownames(Onset_results) <- Scores_to_test

write.csv(Onset_results, file="Results/SZ_ASRB_binary_onset_age_PES_PRS_results.csv", quote = F, row.names = T)

