#############################

## Correlation of PES/PRS with self-reported MHQ phenotypes

## SZ and BIP excluded + controls in the discovery analyses

## Sensitivity analyses

## William Reay (2021)

############################


setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_self_reported_mental_illness/")

## Read in PES/PRS

PES_PRS <- fread("../UKBB_MHQ_no_self_reported_biochem/Entire_EUR_genotyped_UKBB_scores/FINAL_SCORES_PES_PRS/All_SZ_BIP_PES_PRS_all_UKBB.txt", header = T)

PES_PRS <- PES_PRS[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28)]

## Load in MHQ self-report data

MHQ_self_report <- fread("UKBB_self_report_mental_illness_df_SZ_BIP_excl.txt", header = T)

MHQ_self_report <- rename(MHQ_self_report, "IID"="f.eid")

## Merge

Merged_MHQ_scoring <- merge(PES_PRS, MHQ_self_report, by = "IID")

## Define list of scores

Scores_to_test <- as.list(colnames(Merged_MHQ_scoring[,2:15]))