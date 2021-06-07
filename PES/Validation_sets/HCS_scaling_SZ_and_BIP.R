########################################

## Evaluating the proportion of SZ and BIP individuals with elevated PES/PRS relative to the mean of controls in the HCS

## William Reay (2021)

#####################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/")

## Load dependencies

library(dplyr)
library(data.table)
library(easyGgplot2)
library(ggpubr)

## Load HCS data and identify individuals with non-missing K10 psychological distress, CESD scale, and self-reported depression/anxiety

HCS_pheno <- read.csv("../Hunter_cohort/Phenotype_data/mcv1.csv", header = T)

HCS_pheno <- HCS_pheno %>% filter(!is.na(k10w1) & !is.na(CESDw1) & !is.na(Depression_Anxietyw1))

## Remove individuals who meet the following criteria to yield controls for scaling
## i) No self-reported anxiety or depression
## ii) CESD < 20 - PMID: 27182821
## iii) K10 < 17 - PMID: 29698459

HCS_non_psych <- HCS_pheno %>% filter(Depression_Anxietyw1 == 0 &
                                      k10w1 < 17 & CESDw1 < 20)
                                        
## Load SZ scores and append SZ prefix to them

SZ_scores <- fread("Best_DA_gene_PES/HCS_cohort_profiling/SZ_HCS_scores/Combined_SZ_HCS_scores.txt", header = T)

Colnames_all_SZ_scores <- make.unique(names(SZ_scores))

colnames(SZ_scores) <- Colnames_all_SZ_scores

SZ_scores <- SZ_scores %>% rename_all( ~ paste0("SZ_", .x))
SZ_scores <- rename(SZ_scores, "IID"="SZ_IID")

## Load BIP scores and append BIP prefix to them

BIP_scores <- fread("Best_DA_gene_PES/HCS_cohort_profiling/BIP_HCS_scores/Combined_BIP_HCS_scores.txt", header = T)

Colnames_all_BIP_scores <- make.unique(names(BIP_scores))

colnames(BIP_scores) <- Colnames_all_BIP_scores

BIP_scores <- BIP_scores %>% rename_all( ~ paste0("BIP_", .x))
BIP_scores <- rename(BIP_scores, "IID"="BIP_IID")

## Merge all three

HCS_non_psych <- rename(HCS_non_psych, "IID"="electoralId")

SZ_and_BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(SZ_scores, BIP_scores, HCS_non_psych))

## Take the average of the ASRB removed scores in this cohort

Mean_ASRB_vec <- as.vector(SZ_and_BIP_merged[, c(10, 20, 30, 40, 50, 60, 70)] %>% summarise_if(is.numeric, mean))


## SD of the ASRB removed scores in this cohort

SD_ASRB_vec <- as.vector(SZ_and_BIP_merged[, c(10, 20, 30, 40, 50, 60, 70)] %>% summarise_if(is.numeric, sd))

## Take the average of the BIP scores in this cohort

Mean_BIP_vec <- as.vector(SZ_and_BIP_merged[, c(74, 79, 84, 89, 94, 99, 104)] %>% summarise_if(is.numeric, mean))

SD_BIP_vec <- as.vector(SZ_and_BIP_merged[, c(74, 79, 84, 89, 94, 99, 104)] %>% summarise_if(is.numeric, sd))

## Load in ASRB cohort for SZ scaling

All_SZ_scores <- fread("Best_DA_gene_PES/SZ_best_con_or_lib/ASRB_SZ_VALIDATION/Combined_ASRB_replication.txt", header = T, sep = "\t")

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores

## Read in phenotype and covariate data

SZ_pheno <- read.csv("Best_DA_gene_PES/SZ_best_con_or_lib/ASRB_SZ_VALIDATION/Merged_ASRB_case_control.csv", header = T)

## Merge all three df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                    list(All_SZ_scores, SZ_pheno))

## Define only SZ cases

SZ_merged <- SZ_merged %>% filter(PHENOTYPE == 1)


Scores_to_test <- as.list(colnames(SZ_merged[,c(6, 12, 18, 24, 30, 36, 42)]))

## Define function to scale scores to control SD units based on mean and SD in HCS

SD_HC_scale <- function(score_input, control_mean, control_SD, df){
  df[[score_input]] <- (df[[score_input]] - control_mean)/control_SD
  return(df)
}

for (i in 1:7) {
  SD_HC_scale(Scores_to_test[[i]], Mean_ASRB_vec[[i]], SD_ASRB_vec[[i]], SZ_merged)
}


TEST <- mapply(SD_HC_scale, Scores_to_test, Mean_ASRB_vec, SD_ASRB_vec, SZ_merged)


TEST2 <- SD_HC_scale("CACNA1C_PES_AVG", -0.01647076, 0.001225438, SZ_merged)

SZ_merged$Scaled_PRS <- (SZ_merged$SZ_PRS_AVG - Mean_ASRB_vec[[6]])/SD_ASRB_vec[[6]]
SZ_merged$Scaled_CACNA1C_PES <- (SZ_merged$CACNA1C_PES_AVG - Mean_ASRB_vec[[1]])/SD_ASRB_vec[[1]]
SZ_merged$Scaled_FADS1_PES <- (SZ_merged$FADS1_PES_AVG - Mean_ASRB_vec[[2]])/SD_ASRB_vec[[2]]
SZ_merged$Scaled_FES_PES <- (SZ_merged$FES_PES_AVG - Mean_ASRB_vec[[3]])/SD_ASRB_vec[[3]]   
SZ_merged$Scaled_GRIN2A_PES <- (SZ_merged$GRIN2A_PES_AVG - Mean_ASRB_vec[[4]])/SD_ASRB_vec[[4]]
SZ_merged$Scaled_PCCB_PES <- (SZ_merged$PCCB_PES_AVG - Mean_ASRB_vec[[5]])/SD_ASRB_vec[[5]]   
SZ_merged$Scaled_RPS17_PES <- (SZ_merged$RPS17_PES_AVG - Mean_ASRB_vec[[7]])/SD_ASRB_vec[[7]]  

Melted_ASRB_1 <- melt(SZ_merged[,51:53])

Melted_ASRB_2 <- melt(SZ_merged[,54:56])

ggplot(Melted_ASRB_1,aes(x=value,fill=variable))+geom_density()+facet_wrap(variable~., strip.position = "top", nrow=3)+
  xlab("PES (control SD units)") +
  theme_bw() +
  theme(legend.position = "none")

ggplot(Melted_ASRB_2,aes(x=value,fill=variable))+geom_density()+facet_wrap(variable~., strip.position = "top", nrow=3)+
  xlab("PES (control SD units)") +
  theme_bw() +
  theme(legend.position = "none") 

## Plot CACNA1C in the ASRB and HCS

HCS_CACNA1C <- ggplot(SZ_and_BIP_merged,aes(x=SZ_FADS1_PES_NO_ASRB_AVG))+geom_density(fill="grey")+
  xlab("PES (raw average score)") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("HCS SZ CACNA1C PES distribution") +
  xlim(-0.0118837, 0.0152063)

ASRB_CACNA1C <- ggplot(SZ_merged,aes(x=CACNA1C_PES_AVG))+geom_density(fill="#2E9FDF")+
  xlab("PES (raw average score)") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("ASRB SZ CACNA1C PES distribution") +
  xlim(-0.0118837, 0.0152063)

CACNA1C_fig <- ggarrange(HCS_CACNA1C, ASRB_CACNA1C, nrow=2)

## Plot GRIN2A in the ASRB and HCS

HCS_GRIN2A <- ggplot(SZ_and_BIP_merged,aes(x=SZ_GRIN2A_PES_NO_ASRB_AVG))+geom_density(fill="grey")+
  xlab("PES (raw average score)") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("HCS SZ GRIN2A PES distribution") +
  xlim(-2.69052e-06, 6.64431e-06)

ASRB_GRIN2A <- ggplot(SZ_merged,aes(x=GRIN2A_PES_AVG))+geom_density(fill="#E7B800")+
  xlab("PES (raw average score)") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("ASRB SZ GRIN2A PES distribution") +
  xlim(-2.69052e-06, 6.64431e-06) 

GRIN2A_fig <- ggarrange(HCS_GRIN2A, ASRB_GRIN2A, nrow=2)

## Bipolar scaling

All_BIP_scores <- fread("Best_DA_gene_PES/BIP_best_con_or_lib/UKBB_BIP_VALIDATION/Combined_UKBB_BIP_replication_scores.txt", header = T, sep = "\t")

Colnames_all_BIP_scores <- make.unique(names(All_BIP_scores))

colnames(All_BIP_scores) <- Colnames_all_BIP_scores

## Read in phenotype and covariate data

BIP_pheno <- fread("UKBB_phenotypes/BIP_phenotypes/BIP_TEST_UKBB_pheno.txt", header = T)

BIP_cov <- fread("UKBB_phenotypes/BIP_phenotypes/Covariates_TEST_BIP_UKBB.txt", header = T)

## Merge all four df

BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                     list(All_BIP_scores, BIP_pheno, BIP_cov))

BIP_merged$Batch <- as.factor(BIP_merged$Batch)


BIP_merged$Scaled_PRS <- (BIP_merged$BIP_PRS_AVG - Mean_ASRB_vec[[6]])/SD_ASRB_vec[[6]]
BIP_merged$Scaled_CACNA1C_PES <- (BIP_merged$CACNA1C_PES_AVG - Mean_ASRB_vec[[1]])/SD_ASRB_vec[[1]]
BIP_merged$Scaled_FADS1_PES <- (BIP_merged$FADS1_PES_AVG - Mean_ASRB_vec[[2]])/SD_ASRB_vec[[2]]
BIP_merged$Scaled_FES_PES <- (BIP_merged$FES_PES_AVG - Mean_ASRB_vec[[3]])/SD_ASRB_vec[[3]]   
BIP_merged$Scaled_GRIN2A_PES <- (BIP_merged$GRIN2A_PES_AVG - Mean_ASRB_vec[[4]])/SD_ASRB_vec[[4]]
BIP_merged$Scaled_PCCB_PES <- (BIP_merged$PCCB_PES_AVG - Mean_ASRB_vec[[5]])/SD_ASRB_vec[[5]]   
BIP_merged$Scaled_RPS17_PES <- (BIP_merged$RPS17_PES_AVG - Mean_ASRB_vec[[7]])/SD_ASRB_vec[[7]] 


Melted_BIP_1 <- melt(BIP_merged[,61:63])

Melted_BIP_2 <- melt(BIP_merged[,64:65])

ggplot(Melted_BIP_1,aes(x=value,fill=variable))+geom_density()+facet_wrap(variable~., strip.position = "top", nrow=3)+
  xlab("PES (control SD units)") +
  theme_bw() +
  theme(legend.position = "none")

ggplot(Melted_BIP_2,aes(x=value,fill=variable))+geom_density()+facet_wrap(variable~., strip.position = "top", nrow=3)+
  xlab("PES (control SD units)") +
  theme_bw() +
  theme(legend.position = "none") 



