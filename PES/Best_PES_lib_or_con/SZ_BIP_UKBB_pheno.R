#############################

## Defining schizophrenia and bipolar disorder phenotypes in the UKBB

## William Reay (2021)

############################
library(dplyr)
library(data.table)
library(ggplot2)
library(caret)

## Load in ICD-10 and self-reported at baseline

bd <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/Pneumonia_phenotypes/Pneumonia_phenotypes.tab", header = T)

## Read in individuals who have QC EUR genotypes (available for genetic scoring)

GWAS_ind <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/Covariates_pneumonia_UKBB_white_British.tab.txt", header = T)
GWAS_ind <- rename(GWAS_ind, "f.eid"="IID")

## Merge 

Merged_genotype_avail <- merge(bd, GWAS_ind, by="f.eid")

## SELF-REPORT AT BASELINE ##

Self_report_cols <- as.list(paste("f.20002.",0:3, ".", 1:33,sep=""))

select_eid <- Merged_genotype_avail %>%
  select(starts_with("f.eid"))

select_20002 <-  Merged_genotype_avail %>%
  select(starts_with("f.20002."))

Self_reported_df <- bind_cols(select_eid, select_20002)

## Define function to select-self report codes at baseline

Self_report_func <- function(df, code, name) {
  Self_reported <- which(df==code, arr.ind=T)
  Self_reported_row <- Self_reported[,1]
  df$phenotype <- 0
  df[Self_reported_row, ]$phenotype <- 1
  df$phenotype <- as.factor(df$phenotype)
  df <- df %>% select(f.eid, phenotype)
  names(df)[names(df) == "phenotype"] <- name
  return(df)
}

## Extract BIP and SZ data

SZ_baseline_self_report <- Self_report_func(Self_reported_df, "1289", "Self_report_baseline_SZ")

BIP_baseline_self_report <- Self_report_func(Self_reported_df, "1291", "Self_report_baseline_BIP")

Merged_baseline <- merge(SZ_baseline_self_report, BIP_baseline_self_report, by="f.eid")

### ICD-10 PRIMARY OR SECONDARY CODES ###

select_ICD10_primary <- Merged_genotype_avail %>% select(starts_with("f.41202"))

select_ICD10_secondary <- Merged_genotype_avail %>% select(starts_with("f.41204"))

ICD_df <- cbind(select_eid, select_ICD10_primary, select_ICD10_secondary)

## Define list of ICD codes to extract for each disorder

SZ_codes <- as.list(c("F200", "F201", "F202", "F203", "F204",
                      "F205", "F206", "F208", "F209"))

## Extract SZ codes

SZ_ICD_10_select <- which(ICD_df=="F200" | ICD_df=="F201" | ICD_df == "F202" |
                            ICD_df == "F203" | ICD_df == "F204" | ICD_df == "F205" |
                            ICD_df == "F206" | ICD_df == "F208" | ICD_df == "F209", arr.ind = T)

SZ_ICD_10_select_row <- SZ_ICD_10_select[,1]

SZ_ICD_10_select_row <- unique(SZ_ICD_10_select_row)

ICD_df$SZ_ICD10 <- 0
ICD_df[SZ_ICD_10_select_row, ]$SZ_ICD10 <- 1
ICD_df$SZ_ICD10 <- as.factor(ICD_df$SZ_ICD10)
 

## Define list of BIP ICD-19 codes

BIP_codes <- as.list(c("F310", "F311", "F312", "F313", "F314",
                       "F315", "F316", "F317", "F318", "F319") )

BIP_ICD10_select <- which(ICD_df == "F310" | ICD_df == "F311" | ICD_df == "F312" |
                            ICD_df == "F313" | ICD_df == "F314" | ICD_df == "F315" |
                            ICD_df == "F316" | ICD_df == "F317" | ICD_df == "F318" |
                            ICD_df == "F319", arr.ind = T)

BIP_ICD10_select_row <- BIP_ICD10_select[,1]
BIP_ICD10_select_row <- unique(BIP_ICD10_select_row)

ICD_df$BIP_ICD10 <- 0
ICD_df[BIP_ICD10_select_row, ]$BIP_ICD10 <- 1
ICD_df$BIP_ICD10 <- as.factor(ICD_df$BIP_ICD10)

ICD_df <- ICD_df %>% select(f.eid, SZ_ICD10, BIP_ICD10)


## MHQ - SELF REPORT QUESTIONAIRE ##

MHQ_raw <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/MHQ_ever_diagnosed_UKBB.txt", header = T)

## Define individuals who completed the questionaire

MHQ_completed <- MHQ_raw %>% filter(!is.na(f.20400.0.0))

## Merge with available IDs

GWAS_ID <- GWAS_ind %>% select(f.eid)

MHQ_completed <- merge(MHQ_completed, GWAS_ID, by = "f.eid")

## Extract SZ (2) and BIP (10) codes

SZ_MHQ_select <- which(MHQ_completed == 2, arr.ind = T)
SZ_MHQ_select_rows <- SZ_MHQ_select[,1]

MHQ_completed$SZ_MHQ_self_report <- 0
MHQ_completed[SZ_MHQ_select_rows, ]$SZ_MHQ_self_report <- 1
MHQ_completed$SZ_MHQ_self_report <- as.factor(MHQ_completed$SZ_MHQ_self_report)


BIP_MHQ_select <- which(MHQ_completed == 10, arr.ind = T)
BIP_MHQ_select_rows <- BIP_MHQ_select[,1]

MHQ_completed$BIP_MHQ_self_report <- 0
MHQ_completed[BIP_MHQ_select_rows, ]$BIP_MHQ_self_report <- 1
MHQ_completed$BIP_MHQ_self_report <- as.factor(MHQ_completed$BIP_MHQ_self_report)

## Look at the union of all three phenotype types and make list of cases

## SZ ##

Select_MHQ_cases_SZ <- MHQ_completed %>% select(f.eid, SZ_MHQ_self_report, BIP_MHQ_self_report) %>% 
filter(SZ_MHQ_self_report == "1")

Select_ICD_cases_SZ <- ICD_df %>% filter(SZ_ICD10 == "1")

Select_baseline_report_cases_SZ <- Merged_baseline %>% filter(Self_report_baseline_SZ == "1")

## BIP ##

Select_MHQ_cases_BIP <- MHQ_completed %>% select(f.eid, SZ_MHQ_self_report, BIP_MHQ_self_report) %>% 
  filter(BIP_MHQ_self_report == "1")

Select_ICD_cases_BIP <- ICD_df %>% filter(BIP_ICD10 == "1")

Select_baseline_report_cases_BIP <- Merged_baseline %>% filter(Self_report_baseline_BIP == "1")

TE_merge <- merge(Select_ICD_cases_BIP, Select_baseline_report_cases_BIP, by = "f.eid")


## Merge all together

## SZ

Select_MHQ_cases_SZ <- rename(Select_MHQ_cases_SZ, "SZ"="SZ_MHQ_self_report")
Select_MHQ_cases_SZ <- Select_MHQ_cases_SZ %>% select(f.eid, SZ)
Select_baseline_report_cases_SZ <- rename(Select_baseline_report_cases_SZ, "SZ"="Self_report_baseline_SZ")
Select_baseline_report_cases_SZ <- Select_baseline_report_cases_SZ %>% select(f.eid, SZ)
Select_ICD_cases_SZ <- rename(Select_ICD_cases_SZ, "SZ"="SZ_ICD10")
Select_ICD_cases_SZ <- Select_ICD_cases_SZ %>% select(f.eid, SZ)

All_SZ_cases_EUR_unrel_geno_UKBB <- rbind(Select_MHQ_cases_SZ, 
                                          Select_baseline_report_cases_SZ,
                                          Select_ICD_cases_SZ)

All_SZ_cases_EUR_unrel_geno_UKBB <- unique(All_SZ_cases_EUR_unrel_geno_UKBB)

write.table(All_SZ_cases_EUR_unrel_geno_UKBB, file="~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/SZ_phenotypes/All_SZ_unrelated_EUR_genotype_avail_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## BIP

Select_MHQ_cases_BIP <- rename(Select_MHQ_cases_BIP, "BIP"="BIP_MHQ_self_report")
Select_MHQ_cases_BIP <- Select_MHQ_cases_BIP %>% select(f.eid, BIP)
Select_baseline_report_cases_BIP <- rename(Select_baseline_report_cases_BIP, "BIP"="Self_report_baseline_BIP")
Select_baseline_report_cases_BIP <- Select_baseline_report_cases_BIP %>% select(f.eid, BIP)
Select_ICD_cases_BIP <- rename(Select_ICD_cases_BIP, "BIP"="BIP_ICD10")
Select_ICD_cases_BIP <- Select_ICD_cases_BIP %>% select(f.eid, BIP)

All_BIP_cases_EUR_unrel_geno_UKBB <- rbind(Select_MHQ_cases_BIP, 
                                          Select_baseline_report_cases_BIP,
                                          Select_ICD_cases_BIP)

All_BIP_cases_EUR_unrel_geno_UKBB <- unique(All_BIP_cases_EUR_unrel_geno_UKBB)

Overlap <- merge(All_BIP_cases_EUR_unrel_geno_UKBB, All_SZ_cases_EUR_unrel_geno_UKBB, by = "f.eid")

## Define 70% of BIP cases as training set and 30% as test set

trainIndex <- createDataPartition(All_BIP_cases_EUR_unrel_geno_UKBB$f.eid,p=0.7,list=FALSE)


train_BIP <- All_BIP_cases_EUR_unrel_geno_UKBB[trainIndex, ]
test_BIP <- All_BIP_cases_EUR_unrel_geno_UKBB[-trainIndex, ]

## Ensure no overlap

Merge_train_test <- merge(train_BIP, test_BIP, by = "f.eid")

write.table(train_BIP, file="~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/BIP_phenotypes/TRAINING_BIP_unrelated_EUR_genotype_avail_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

write.table(test_BIP, file="~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/BIP_phenotypes/TEST_BIP_unrelated_EUR_genotype_avail_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## Define controls without relevant codes and no self-reported conditions in the MHQ

No_MH_select <- which(!is.na(MHQ_completed[, 3:17]), arr.ind = T)
No_MH_select_rows <- No_MH_select[,1]
No_MH_select_rows <- unique(No_MH_select_rows)

MHQ_completed$No_mental_health <- 0
MHQ_completed[No_MH_select_rows, ]$No_mental_health <- 1
MHQ_completed$No_mental_health <- as.factor(MHQ_completed$No_mental_health)

No_report_MH <- MHQ_completed %>% select(f.eid, No_mental_health) %>% filter(No_mental_health == 0)

## Merge with SZ and BIP baseline characteristics to check for no overlap

Merged_SZ_no_MH <- merge(No_report_MH, All_SZ_cases_EUR_unrel_geno_UKBB,
                         by = "f.eid")

Merged_BIP_no_MH <- merge(No_report_MH, All_BIP_cases_EUR_unrel_geno_UKBB,
                         by = "f.eid")


## Remove the three overlapping individuals

No_report_MH_cases_removed <- No_report_MH %>% filter(f.eid != 2839592 &
                                                        f.eid != 3139224 &
                                                        f.eid != 4964352)

No_report_MH_cases_removed$SZ <- 0
No_report_MH_cases_removed$BIP <- 0

## Select double the number of controls for each cohort
## SZ - N = 1262
## BIP training - N = 2318
## BIP test - N = 996

SZ_controls <- No_report_MH_cases_removed[1:1262, ]

SZ_controls <- SZ_controls %>% select(f.eid, SZ)

BIP_training_controls <- No_report_MH_cases_removed[25000:27317, ]

BIP_training_controls <- BIP_training_controls %>% select(f.eid, BIP)

BIP_test_controls <- No_report_MH_cases_removed[40000:40995, ]

BIP_test_controls <- BIP_test_controls %>% select(f.eid, BIP)

## Bind with relevant case df and combine with covariates

SZ_UKBB_full <- rbind(All_SZ_cases_EUR_unrel_geno_UKBB, SZ_controls)

SZ_UKBB_full <- merge(SZ_UKBB_full, GWAS_ind, by = "f.eid")

write.table(SZ_UKBB_full, file="~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/SZ_phenotypes/SZ_case_double_N_random_control_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

BIP_UKBB_training <- rbind(train_BIP, BIP_training_controls)

BIP_training_UKBB_full <- merge(BIP_UKBB_training, GWAS_ind, by = "f.eid")

write.table(BIP_training_UKBB_full, file="~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/BIP_phenotypes/BIP_TRAINING_double_N_random_control_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

BIP_UKBB_test <- rbind(test_BIP, BIP_test_controls)

BIP_test_UKBB_full <- merge(BIP_UKBB_test, GWAS_ind, by = "f.eid")

write.table(BIP_test_UKBB_full, file="~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/BIP_phenotypes/BIP_TEST_double_N_random_control_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## Test for overlap

SZ <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/SZ_phenotypes/SZ_case_double_N_random_control_UKBB.txt", header = T)

BIP_train <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/BIP_phenotypes/BIP_TRAINING_double_N_random_control_UKBB.txt", header = T)

BIP_test <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/BIP_phenotypes/BIP_TEST_double_N_random_control_UKBB.txt", header = T)

Final_ov_test <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "f.eid"),
                                         list(SZ, BIP_test,
                                              BIP_train))


## Demographic characteristics

SZ_demo <- glm(SZ ~ Sex + Age, family="binomial", data = SZ_UKBB_full)

BIP_train_demo <- glm(BIP ~ Sex + Age, family="binomial", data = BIP_training_UKBB_full)

BIP_test_demo <- glm(BIP ~ Sex + Age, family="binomial", data = BIP_test_UKBB_full)
