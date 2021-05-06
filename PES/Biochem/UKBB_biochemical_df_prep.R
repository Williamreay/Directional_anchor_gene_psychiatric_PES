#############################

## Defining the cohort for the biochemical pheWAS

## MHQ participants without a major self-reported mental illness and also not included in any of the SZ/BIP cohorts

## William Reay (2021)

############################

library(dplyr)
library(data.table)

## Read in individuals who have QC EUR genotypes (available for genetic scoring)

GWAS_ind <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/Covariates_pneumonia_UKBB_white_British.tab.txt", header = T)
GWAS_ind <- rename(GWAS_ind, "f.eid"="IID")

## MHQ - SELF REPORT QUESTIONAIRE ##

MHQ_raw <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/MHQ_ever_diagnosed_UKBB.txt", header = T)

## Define individuals who completed the questionaire

MHQ_completed <- MHQ_raw %>% filter(!is.na(f.20400.0.0))

## Define controls without relevant codes and no self-reported conditions in the MHQ

No_MH_select <- which(!is.na(MHQ_completed[, 3:17]), arr.ind = T)
No_MH_select_rows <- No_MH_select[,1]
No_MH_select_rows <- unique(No_MH_select_rows)

MHQ_completed$No_mental_health <- 0
MHQ_completed[No_MH_select_rows, ]$No_mental_health <- 1
MHQ_completed$No_mental_health <- as.factor(MHQ_completed$No_mental_health)

No_report_MH <- MHQ_completed %>% select(f.eid, No_mental_health) %>% filter(No_mental_health == 0)

## Merge with available GWAS IDs

No_MH_MHQ_completed <- merge(No_report_MH, GWAS_ind, by = "f.eid")

## Remove any individual who was in the SZ and/or BIP cohorts

Overlap <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/BIP_and_SZ_IDs.txt", header = T)

Overlap <- rename(Overlap, "f.eid"="IID", "FID"="#FID")

Biochem_cohort_merged <- merge(No_MH_MHQ_completed, Overlap, by = "f.eid", all.x = T)

Final_Biochem_cohort_merged <- Biochem_cohort_merged %>% filter(is.na(FID))

## Read in biochemical df

Raw_biochem <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/UKBB_biochemical_pheWAS_IDs/Meausred_biochemical_data_UKBB.txt", header = T)

## Merge

Final_Biochem_cohort_merged <- merge(Final_Biochem_cohort_merged, Raw_biochem, by = "f.eid")

write.table(Final_Biochem_cohort_merged, file="~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/UKBB_biochemical_pheWAS_IDs/Measured_MHQ_controls_biochem.txt",
            sep = "\t", row.names = F, quote = F)
