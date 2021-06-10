#############################

## Defining remaining self-reported mental health phenotypes in the MHQ

## SZ and BIP excluded + controls in the discovery analyses

## William Reay (2021)

############################

library(dplyr)
library(data.table)

## MHQ - SELF REPORT QUESTIONAIRE ##

MHQ_raw <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/MHQ_ever_diagnosed_UKBB.txt", header = T)

## Define individuals who completed the questionaire

MHQ_completed <- MHQ_raw %>% filter(!is.na(f.20400.0.0))

## Merge with individuals with genotype data available

GWAS_ind <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/Covariates_pneumonia_UKBB_white_British.tab.txt", header = T)
GWAS_ind <- rename(GWAS_ind, "f.eid"="IID")

MHQ_completed <- merge(MHQ_completed, GWAS_ind, by="f.eid")

## Remove any individual who was in the SZ and/or BIP cohorts

Overlap <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/BIP_and_SZ_IDs.txt", header = T)

Overlap <- rename(Overlap, "f.eid"="IID", "FID"="#FID")

MHQ_completed_no_SZ_BIP <- merge(MHQ_completed, Overlap, by = "f.eid", all.x = T)

MHQ_completed_no_SZ_BIP<- MHQ_completed_no_SZ_BIP %>% filter(is.na(FID))


## Define control cohort

## Define controls without relevant codes and no self-reported conditions in the MHQ

No_MH_select <- which(!is.na(MHQ_completed_no_SZ_BIP[, 3:17]), arr.ind = T)
No_MH_select_rows <- No_MH_select[,1]
No_MH_select_rows <- unique(No_MH_select_rows)

MHQ_completed_no_SZ_BIP$No_mental_health <- 0
MHQ_completed_no_SZ_BIP[No_MH_select_rows, ]$No_mental_health <- 1
MHQ_completed_no_SZ_BIP$No_mental_health <- as.factor(MHQ_completed_no_SZ_BIP$No_mental_health)

## Retain only the self_reported rows

select_eid <- MHQ_completed_no_SZ_BIP %>%
  select(starts_with("f.eid"))

select_20544 <-  MHQ_completed_no_SZ_BIP %>%
  select(starts_with("f.20544."))

Self_reported_df <- bind_cols(select_eid, select_20544)

## Pass function to define phenotypic categories in controls such that it can merged with the 'case df'

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


## List of codes

Codes <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_self_reported_mental_illness/Self_report_MH_coding.txt", header = T)

Codes_input <- as.list(Codes$coding)

Name_input <- as.list(Codes$Input_meaning)

MHQ_completed_no_SZ_BIP_cols_added <- mapply(Self_report_func, code = Codes_input,
                       name = Name_input, MoreArgs = list(Self_reported_df),
                       SIMPLIFY = FALSE)

MHQ_completed_no_SZ_BIP_cols_added <- as.data.frame(MHQ_completed_no_SZ_BIP_cols_added)

## Final merge

MHQ_completed_no_SZ_BIP_cols_added <- merge(MHQ_completed_no_SZ_BIP, MHQ_completed_no_SZ_BIP_cols_added, by = "f.eid")


## Write ouput

write.table(MHQ_completed_no_SZ_BIP_cols_added, 
            file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_self_reported_mental_illness/UKBB_self_report_mental_illness_df_SZ_BIP_excl.txt",
            sep = "\t", row.names = F, quote = F)



