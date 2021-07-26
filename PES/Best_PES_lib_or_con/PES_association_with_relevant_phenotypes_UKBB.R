################################

## Identifying PES relevant phenotypes for UKBB participants who did not self-report mental illness in the MHQ

## SZ and BIP cohorts excluded

## William Reay (2021)

#################################

library(dplyr)
library(data.table)

## Read in the MHQ no-self reported cohort (SZ and BIP control cohorts excluded) along with the scores

MHQ_cohort <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem/Merged_all_scores_all_biochem.txt", header = T)

## Merge with self-reported and ICD-10 df

bd <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/Pneumonia_phenotypes/Pneumonia_phenotypes.tab", header = T)

bd <- rename(bd, "IID"="f.eid")

MHQ_cohort_merged <- merge(MHQ_cohort, bd, by ="IID")

## Define the following phenotypes for each of the PES using self-report and/or ICD-10 data

## CACNA1C PES - atrial fibrilation - self report (1471, 1483), ICD-10 (I480, I481, I482, I483, I484, I489)
## FADS1 PES - Hypercholesteremia - self-report (1473), do not include familal ICD-10 codes
## FES PES - chronic immune thrombocytopenia (ITP) - ICD-10 (D693)
## GRIN2A PES - COPD - self-report (1112) and ICD-10 (J440, J441, J448, J449)
## PCCB PES - vitamin B deficency - ICD-10 (E359) and (D510, D511, D512, D513, D518, D519)
## RPS17 PES - malaria, however, will not include

## SELF-REPORT AT BASELINE ##

Self_report_cols <- as.list(paste("f.20002.",0:3, ".", 1:33,sep=""))

select_eid <- MHQ_cohort_merged %>%
  select(starts_with("IID"))

select_20002 <-  MHQ_cohort_merged %>%
  select(starts_with("f.20002."))

Self_reported_df <- bind_cols(select_eid, select_20002)

## Define function to select-self report codes

Self_report_func <- function(df, code, name) {
  Self_reported <- which(df==code, arr.ind=T)
  Self_reported_row <- Self_reported[,1]
  df$phenotype <- 0
  df[Self_reported_row, ]$phenotype <- 1
  df$phenotype <- as.factor(df$phenotype)
  df <- df %>% select(IID, phenotype)
  names(df)[names(df) == "phenotype"] <- name
  return(df)
}

## Extract dataframe for each relevant phenotype

# Atrial fibrillation
AFib_self_report <- Self_report_func(Self_reported_df, "1471", "AFib")
AFlut_self_report <- Self_report_func(Self_reported_df, "1483", "AFlut")

AF <- merge(AFib_self_report, AFlut_self_report, by="IID")
AF$AF_full <- ifelse(AF$AFib == 1 | AF$AFlut == 1, 1, 0)

# Hypercholesteremia
High_chol <- Self_report_func(Self_reported_df, "1473", "High_chol")

# COPD
COPD <- Self_report_func(Self_reported_df, "1112", "COPD")

## Merge all self-report

Self_report_MHQ_incl <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                               list(MHQ_cohort_merged, AF, High_chol, COPD))

## ICD10 inpatient records ##

select_ICD10_primary <- MHQ_cohort_merged %>% select(starts_with("f.41202"))

select_ICD10_secondary <- MHQ_cohort_merged %>% select(starts_with("f.41204"))

ICD_df <- cbind(select_eid, select_ICD10_primary, select_ICD10_secondary)

## Extract Atrial fibrilation codes: ICD-10 (I480, I481, I482, I483, I484, I489)

AF_ICD_10_select <- which(ICD_df=="I480" | ICD_df=="I481" | ICD_df == "I482" |
                            ICD_df == "I483" | ICD_df == "I484" | ICD_df == "I489", arr.ind = T)

AF_ICD_10_select_row <- AF_ICD_10_select[,1]

AF_ICD_10_select_row <- unique(AF_ICD_10_select_row)

ICD_df$AF_ICD10 <- 0
ICD_df[AF_ICD_10_select_row, ]$AF_ICD10 <- 1
ICD_df$AF_ICD10 <- as.factor(ICD_df$AF_ICD10)

## Extract chronic immune thrombocytopenia (ITP) code D693

ITP_ICD_10_select <- which(ICD_df=="D693", arr.ind = T)

ITP_ICD_10_select_row <- ITP_ICD_10_select[,1]

ITP_ICD_10_select_row <- unique(ITP_ICD_10_select_row)

ICD_df$ITP_ICD10 <- 0
ICD_df[ITP_ICD_10_select_row, ]$ITP_ICD10 <- 1
ICD_df$ITP_ICD10 <- as.factor(ICD_df$ITP_ICD10)

## Extract COPD - J440, J441, J448, J449

COPD_ICD_10_select <- which(ICD_df=="J440" | ICD_df=="J441" | ICD_df == "J448" |
                            ICD_df == "J449", arr.ind = T)

COPD_ICD_10_select_row <- COPD_ICD_10_select[,1]

COPD_ICD_10_select_row <- unique(COPD_ICD_10_select_row)

ICD_df$COPD_ICD10 <- 0
ICD_df[COPD_ICD_10_select_row, ]$COPD_ICD10 <- 1
ICD_df$COPD_ICD10 <- as.factor(ICD_df$COPD_ICD10)

## Extract vitamin B deficiency related codes - ICD-10 (E539) and (D510, D511, D512, D513, D518, D519)

VitB_ICD_10_select <- which(ICD_df=="E539" | ICD_df=="D510" | ICD_df == "D511" |
                              ICD_df == "D512" | ICD_df == "D513" | ICD_df == "D518" | ICD_df == "D519", arr.ind = T)

VitB_ICD_10_select_row <- VitB_ICD_10_select[,1]

VitB_ICD_10_select_row <- unique(VitB_ICD_10_select_row)

ICD_df$VitB_def_ICD10 <- 0
ICD_df[VitB_ICD_10_select_row, ]$VitB_def_ICD10 <- 1

## Merge ICD10 with rest 

Final_MHQ_df <- merge(Self_report_MHQ_incl, ICD_df, by="IID")

## Make final phenotype definition for AF and COPD

Final_MHQ_df$AF_final <- ifelse(Final_MHQ_df$AF_full == 1 | Final_MHQ_df$AF_ICD10 == 1, 1, 0)
Final_MHQ_df$COPD_final <- ifelse(Final_MHQ_df$COPD == 1 | Final_MHQ_df$COPD_ICD10 == 1, 1, 0)

## Final prevalences
## AF = 9715
## High cholesterol = 8543
## ITP = 64
## COPD = 761
## Vitamin B defiency = 103

## Scale PES/PRS

Final_MHQ_df[,c(478,480,482,484,486,488,490,492,496,494,498,500,502,504)] <- lapply(Final_MHQ_df[,c(478,480,482,484,486,488,490,492,496,494,498,500,502,504)], function(x) c(scale(x)))

Final_MHQ_df$Batch <- as.factor(Final_MHQ_df$Batch)

## Test each PES set with its respective disorder

## CACNA1C PES and AF

SZ_CACNA1C <- glm(AF_final ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 +  SZ_CACNA1C_PES + Batch, family = "binomial", data = Final_MHQ_df)

BIP_CACNA1C <- glm(AF_final ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 +  BIP_CACNA1C_PES + Batch, family = "binomial", data = Final_MHQ_df)

## FADS1 PES and high cholesterol

SZ_FADS1 <- glm(High_chol ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 +  SZ_FADS1_PES + Batch, family = "binomial", data = Final_MHQ_df)

BIP_FADS1 <- glm(High_chol ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                     PC8 + PC9 + PC10 +  BIP_FADS1_PES + Batch, family = "binomial", data = Final_MHQ_df)

## FES PES and ITP

SZ_FES <- glm(ITP_ICD10 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                  PC8 + PC9 + PC10 +  SZ_FES_PES + Batch, family = "binomial", data = Final_MHQ_df)

BIP_FES <- glm(ITP_ICD10 ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                   PC8 + PC9 + PC10 +  BIP_FES_PES + Batch, family = "binomial", data = Final_MHQ_df)

## GRIN2A PES and COPD

SZ_GRIN2A <- glm(COPD_final ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                PC8 + PC9 + PC10 +  SZ_GRIN2A_PES + Batch, family = "binomial", data = Final_MHQ_df)

BIP_GRIN2A <- glm(COPD_final ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                 PC8 + PC9 + PC10 +  BIP_GRIN2A_PES + Batch, family = "binomial", data = Final_MHQ_df)
