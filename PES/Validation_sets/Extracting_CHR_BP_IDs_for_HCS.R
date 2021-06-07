####################################

## Extracting the score info and output scoring files for best performing PES/PRS

## CHR:BP IDs for HCS

## William Reay (2021)

####################################

## SZ PGC3 ALL ##

library(data.table)
library(dplyr)

## SZ ##

## PRS - penalised regression ## 

PRS <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_PRS.validate.rds")

PRS_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_PRS.lassosum.pipeline.rds")

## Summary statistics for PRS

PRS_weights <- cbind(PRS_pipeline$sumstats, PRS$best.beta)

SZ_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/QC_sumstats/SZ_PGC3_ALL_QC_sumstats.txt.gz", header = T)

SZ_sumstats$order <-  seq.int(nrow(SZ_sumstats))

Merged_PRS <- merge(PRS_weights, SZ_sumstats, by = "order")

PRS_output <- Merged_PRS %>% select('CHR:BP', A1.x, 'PRS$best.beta')

PRS_output <- rename(PRS_output, "A1"="A1.x", "SZ_PRS_ALL"="PRS$best.beta", "SNP"="CHR:BP")

write.table(PRS_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_SZ_PRS_ALL_PGC3.txt",
            row.names = F, quote = F)

## CACNA1C - C+T (conservative boundaries) ## 

CACNA1C_SNPs_to_include <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/UKBB_SZ/SZ_UKBB_SZ_PGC3_ALL_CACNA1C_network.snp", header = T)

CACNA1C_SNPs_to_include <- CACNA1C_SNPs_to_include %>% filter(P < 0.005)

CACNA1C_sum_stats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/SZ_PGC3_ALL_CACNA1C_network.txt", header = T)

CACNA1C_sum_stats$Beta <- log(CACNA1C_sum_stats$OR)

Merged_CACNA1C <- merge(CACNA1C_sum_stats, CACNA1C_SNPs_to_include, by = "SNP")

CACNA1C_output <- Merged_CACNA1C %>% select('CHR:BP', A1, Beta)

CACNA1C_output <- rename(CACNA1C_output, "CACNA1C_PES_ALL"="Beta", "SNP"="CHR:BP")

write.table(CACNA1C_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_CACNA1C_PES_SZ_cons_C_T.txt",
            row.names = F, quote = F)

## FADS1 - C+T (conservative boundaries) ## 

FADS1_SNPs_to_include <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/UKBB_SZ/SZ_UKBB_SZ_PGC3_ALL_FADS1_network.snp", header = T)

FADS1_SNPs_to_include <- FADS1_SNPs_to_include %>% filter(P < 0.05)

FADS1_sum_stats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/SZ_PGC3_ALL_FADS1_network.txt", header = T)

FADS1_sum_stats$Beta <- log(FADS1_sum_stats$OR)

Merged_FADS1 <- merge(FADS1_sum_stats, FADS1_SNPs_to_include, by = "SNP")

FADS1_output <- Merged_FADS1 %>% select("CHR:BP", A1, Beta)

FADS1_output <- rename(FADS1_output, "FADS1_PES_ALL"="Beta", "SNP"="CHR:BP")

write.table(FADS1_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_FADS1_PES_SZ_cons_C_T.txt",
            row.names = F, quote = F)

## FES - C+T (liberal) ##

FES_SNPs_to_include <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/SZ_UKBB_liberal/SZ_UKBB_FES.snp", header = T)

FES_SNPs_to_include <- FES_SNPs_to_include %>% filter(P < 0.05)

FES_sum_stats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/Liberal_boundaries_pathway_specific_variants/Liberal_bound_SZ_all_PGC3_FES", header = T)

FES_sum_stats$Beta <- log(FES_sum_stats$OR)

Merged_FES <- merge(FES_sum_stats, FES_SNPs_to_include, by = "SNP")

FES_output <- Merged_FES %>% select('CHR:BP', A1, Beta)

FES_output <- rename(FES_output, "FES_PES_ALL"="Beta", "SNP"="CHR:BP")

write.table(FES_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_FES_PES_SZ_cons_C_T.txt",
            row.names = F, quote = F)

## GRIN2A - penalised regression ##

GRIN2A <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/SZ_UKBB_lassosum_liberal/SZ_UKBB_LASSO_GRIN2A.validate.rds")

GRIN2A_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/SZ_UKBB_lassosum_liberal/SZ_UKBB_LASSO_GRIN2A.lassosum.pipeline.rds")

## Summary statistics for PRS

GRIN2A_weights <- cbind(GRIN2A_pipeline$sumstats, GRIN2A$best.beta)

NZ <- GRIN2A_weights %>% filter(GRIN2A$best.beta != 0)

GRIN2A_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/Liberal_boundaries_pathway_specific_variants/Liberal_bound_SZ_all_PGC3_GRIN2A", header = T)

GRIN2A_sumstats$order <-  seq.int(nrow(GRIN2A_sumstats))

Merged_GRIN2A <- merge(GRIN2A_weights, GRIN2A_sumstats, by = "order")

GRIN2A_output <- Merged_GRIN2A %>% select("CHR:BP", A1.x, 'GRIN2A$best.beta')

GRIN2A_output <- rename(GRIN2A_output, "GRIN2A_PES_ALL"="GRIN2A$best.beta", "A1"="A1.x")

write.table(GRIN2A_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_GRIN2A_PES_penalised_reg_weights.txt",
            row.names = F, quote = F)

## PCCB - penalised regression ##

PCCB <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_PCCB_network.validate.rds")

PCCB_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_PCCB_network.lassosum.pipeline.rds")

## Summary statistics for PRS

PCCB_weights <- cbind(PCCB_pipeline$sumstats, PCCB$best.beta)

PCCB_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/SZ_PGC3_ALL_PCCB_network.txt", header = T)

PCCB_sumstats$order <-  seq.int(nrow(PCCB_sumstats))

Merged_PCCB <- merge(PCCB_weights, PCCB_sumstats, by = "order")

PCCB_output <- Merged_PCCB %>% select("CHR:BP", A1.x, 'PCCB$best.beta')

PCCB_output <- rename(PCCB_output, "PCCB_PES_ALL"="PCCB$best.beta", "A1"="A1.x")

write.table(PCCB_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_PCCB_PES_penalised_reg_weights.txt",
            row.names = F, quote = F)

## RPS17 - penalised regression (conservative) ##

RPS17 <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_RPS17_network.validate.rds")

RPS17_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_RPS17_network.lassosum.pipeline.rds")

## Summary statistics for PRS

RPS17_weights <- cbind(RPS17_pipeline$sumstats, RPS17$best.beta)

RPS17_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/SZ_PGC3_ALL_RPS17_network.txt", header = T)

RPS17_sumstats$order <-  seq.int(nrow(RPS17_sumstats))

Merged_RPS17 <- merge(RPS17_weights, RPS17_sumstats, by = "order")

RPS17_output <- Merged_RPS17 %>% select("CHR:BP", A1.x, 'RPS17$best.beta')

RPS17_output <- rename(RPS17_output, "RPS17_PES_ALL"="RPS17$best.beta", "A1"="A1.x")

write.table(RPS17_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_RPS17_PES_penalised_reg_weights.txt",
            row.names = F, quote = F)

## SZ - ASRB removed ##

## PRS - penalised regression ## 

PRS <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/ASRB_removed_lassosum_SZ_UKBB_TRAIN/SZ_UKBB_LASSO_PRS_ASRB_removed.validate.rds")

PRS_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/ASRB_removed_lassosum_SZ_UKBB_TRAIN/SZ_UKBB_LASSO_PRS_ASRB_removed.lassosum.pipeline.rds")

## Summary statistics for PRS

PRS_weights <- cbind(PRS_pipeline$sumstats, PRS$best.beta)

SZ_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/QC_sumstats/SZ_PGC3_NO_ASRB_QC_sumstats.txt.gz", header = T)

SZ_sumstats$order <-  seq.int(nrow(SZ_sumstats))

Merged_PRS <- merge(PRS_weights, SZ_sumstats, by = "order")

Colnames_PRS <- make.unique(names(Merged_PRS))

colnames(Merged_PRS) <- Colnames_PRS

PRS_output <- Merged_PRS %>% select('CHR:BP', A1.x, 'PRS$best.beta')

PRS_output <- rename(PRS_output, "A1"="A1.x", "SZ_PRS_NO_ASRB"="PRS$best.beta", "SNP"="CHR:BP")

write.table(PRS_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_ASRB_excluded_SZ_PRS_penalised_reg_weights.txt",
            row.names = F, quote = F)

## CACNA1C - C+T (conservative boundaries) ## 

CACNA1C_SNPs_to_include <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/UKBB_SZ/SZ_UKBB_SZ_PGC3_ALL_CACNA1C_network.snp", header = T)

CACNA1C_SNPs_to_include <- CACNA1C_SNPs_to_include %>% filter(P < 0.005)

CACNA1C_sum_stats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/NO_ASRB_SZ_PGC3_CACNA1C", header = T)

CACNA1C_sum_stats$Beta <- log(CACNA1C_sum_stats$OR)

Merged_CACNA1C <- merge(CACNA1C_sum_stats, CACNA1C_SNPs_to_include, by = "SNP")

CACNA1C_output <- Merged_CACNA1C %>% select('CHR:BP', A1, Beta)

CACNA1C_output <- rename(CACNA1C_output, "CACNA1C_PES_NO_ASRB"="Beta", "SNP"="CHR:BP")

write.table(CACNA1C_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_ASRB_excluded_CACNA1C_PES_SZ_cons_C_T.txt",
            row.names = F, quote = F)

## FADS1 - C+T (conservative boundaries) ## 

FADS1_SNPs_to_include <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/UKBB_SZ/SZ_UKBB_SZ_PGC3_ALL_FADS1_network.snp", header = T)

FADS1_SNPs_to_include <- FADS1_SNPs_to_include %>% filter(P < 0.05)

FADS1_sum_stats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/NO_ASRB_SZ_PGC3_FADS1", header = T)

FADS1_sum_stats$Beta <- log(FADS1_sum_stats$OR)

Merged_FADS1 <- merge(FADS1_sum_stats, FADS1_SNPs_to_include, by = "SNP")

FADS1_output <- Merged_FADS1 %>% select('CHR:BP', A1, Beta)

FADS1_output <- rename(FADS1_output, "FADS1_PES_NO_ASRB"="Beta", "SNP"="CHR:BP")

write.table(FADS1_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_ASRB_removed_FADS1_PES_SZ_cons_C_T.txt",
            row.names = F, quote = F)

## FES - C+T (liberal) ##

FES_SNPs_to_include <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/SZ_UKBB_liberal/SZ_UKBB_FES.snp", header = T)

FES_SNPs_to_include <- FES_SNPs_to_include %>% filter(P < 0.05)

FES_sum_stats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/Liberal_boundaries_pathway_specific_variants/Liberal_bound_NO_ASRB_PGC3_FES", header = T)

FES_sum_stats$Beta <- log(FES_sum_stats$OR)

Merged_FES <- merge(FES_sum_stats, FES_SNPs_to_include, by = "SNP")

FES_output <- Merged_FES %>% select("CHR:BP", A1, Beta)

FES_output <- rename(FES_output, "FES_PES_NO_ASRB"="Beta", "SNP"="CHR:BP")

write.table(FES_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_ASRB_excluded_FES_PES_SZ_cons_C_T.txt",
            row.names = F, quote = F)

## GRIN2A - penalised regression ##

GRIN2A <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/ASRB_GRIN2A_lassosum_training/ASRB_removed_SZ_UKBB_LASSO_GRIN2A.validate.rds")

GRIN2A_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/ASRB_GRIN2A_lassosum_training/ASRB_removed_SZ_UKBB_LASSO_GRIN2A.lassosum.pipeline.rds")

## Summary statistics for PRS

GRIN2A_weights <- cbind(GRIN2A_pipeline$sumstats, GRIN2A$best.beta)

NZ <- GRIN2A_weights %>% filter(GRIN2A$best.beta != 0)

GRIN2A_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/Liberal_boundaries_pathway_specific_variants/Liberal_bound_NO_ASRB_PGC3_GRIN2A", header = T)

GRIN2A_sumstats$order <-  seq.int(nrow(GRIN2A_sumstats))

Merged_GRIN2A <- merge(GRIN2A_weights, GRIN2A_sumstats, by = "order")

GRIN2A_output <- Merged_GRIN2A %>% select('CHR:BP', A1.x, 'GRIN2A$best.beta')

GRIN2A_output <- rename(GRIN2A_output, "GRIN2A_PES_NO_ASRB"="GRIN2A$best.beta", "A1"="A1.x", "SNP"="CHR:BP")

write.table(GRIN2A_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_ASRB_excluded_GRIN2A_PES_penalised_reg_weights.txt",
            row.names = F, quote = F)

## PCCB - penalised regression ##

PCCB <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_PCCB_network.validate.rds")

PCCB_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_PCCB_network.lassosum.pipeline.rds")

## Summary statistics for PRS

PCCB_weights <- cbind(PCCB_pipeline$sumstats, PCCB$best.beta)

PCCB_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/NO_ASRB_SZ_PGC3_PCCB", header = T)

PCCB_sumstats$order <-  seq.int(nrow(PCCB_sumstats))

Merged_PCCB <- merge(PCCB_weights, PCCB_sumstats, by = "order")

PCCB_output <- Merged_PCCB %>% select('CHR:BP', A1.x, 'PCCB$best.beta')

PCCB_output <- rename(PCCB_output, "SNP"="CHR:BP", "PCCB_PES_NO_ASRB"="PCCB$best.beta", "A1"="A1.x")

write.table(PCCB_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_ASRB_excluded_PCCB_PES_penalised_reg_weights.txt",
            row.names = F, quote = F)

## RPS17 - penalised regression (conservative) ##

RPS17 <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_RPS17_network.validate.rds")

RPS17_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_RPS17_network.lassosum.pipeline.rds")

## Summary statistics for PRS

RPS17_weights <- cbind(RPS17_pipeline$sumstats, RPS17$best.beta)

RPS17_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/NO_ASRB_SZ_PGC3_RPS17", header = T)

RPS17_sumstats$order <-  seq.int(nrow(RPS17_sumstats))

Merged_RPS17 <- merge(RPS17_weights, RPS17_sumstats, by = "order")

RPS17_output <- Merged_RPS17 %>% select('CHR:BP', A1.x, 'RPS17$best.beta')

RPS17_output <- rename(RPS17_output, "RPS17_PES_NO_ASRB"="RPS17$best.beta", "A1"="A1.x", "SNP"="CHR:BP")

write.table(RPS17_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_ASRB_excluded_RPS17_PES_penalised_reg_weights.txt",
            row.names = F, quote = F)

## BIP ##

## PRS - penalised regression ## 

PRS <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_BIP_UKBB_TRAIN/BIP_TRAIN_UKBB_LASSO_BIP_TRAIN_PRS.validate.rds")

PRS_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_BIP_UKBB_TRAIN/BIP_TRAIN_UKBB_LASSO_BIP_TRAIN_PRS.lassosum.pipeline.rds")

## Summary statistics for PRS

PRS_weights <- cbind(PRS_pipeline$sumstats, PRS$best.beta)

BIP_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/QC_sumstats/BIP_QC_sumstats.txt.gz", header = T)

BIP_sumstats$order <-  seq.int(nrow(BIP_sumstats))

Merged_PRS <- merge(PRS_weights, BIP_sumstats, by = "order")

PRS_output <- Merged_PRS %>% select('CHR:BP', A1.x, 'PRS$best.beta')

PRS_output <- rename(PRS_output, "A1"="A1.x", "BIP_PRS"="PRS$best.beta", "SNP"="CHR:BP")

write.table(PRS_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_Penalised_regression_BIP_gw_PRS.txt",
            row.names = F, quote = F)

## CACNA1C - Penalised regression (liberal) ## 

CACNA1C <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/BIP_TRAIN_lassosum_liberal/BIP_TRAIN_UKBB_LASSO_CACNA1C.validate.rds")

CACNA1C_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/BIP_TRAIN_lassosum_liberal/BIP_TRAIN_UKBB_LASSO_CACNA1C.lassosum.pipeline.rds")

## Summary statistics for CACNA1C

CACNA1C_weights <- cbind(CACNA1C_pipeline$sumstats, CACNA1C$best.beta)

CACNA1C_NZ <- CACNA1C_weights %>% filter(CACNA1C$best.beta != 0)

CACNA1C_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/Liberal_boundaries_pathway_specific_variants/Liberal_bound_BIP_CACNA1C", header = T)

CACNA1C_sumstats$order <-  seq.int(nrow(CACNA1C_sumstats))

Merged_CACNA1C <- merge(CACNA1C_weights, CACNA1C_sumstats, by = "order")

CACNA1C_output <- Merged_CACNA1C %>% select("CHR:BP", A1.x, 'CACNA1C$best.beta')

CACNA1C_output <- rename(CACNA1C_output, "A1"="A1.x", "CACNA1C_PES"="CACNA1C$best.beta", "SNP"="CHR:BP")

write.table(CACNA1C_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_CACNA1C_PES_BIP_penalised_reg_lib.txt",
            row.names = F, quote = F)

## FADS1 - C+T (conservative)

FADS1_SNPs_to_include <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/BIP_TRAIN_UKBB/BIP_TRAIN_UKBB_BIP_TRAIN_FADS1_network.snp", header = T)

FADS1_sum_stats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/BIP_FADS1_network.txt", header = T)

FADS1_sum_stats$Beta <- log(FADS1_sum_stats$OR)

Merged_FADS1 <- merge(FADS1_sum_stats, FADS1_SNPs_to_include, by = "SNP")

FADS1_output <- Merged_FADS1 %>% select('CHR:BP', A1, Beta)

FADS1_output <- rename(FADS1_output, "FADS1_PES"="Beta", "SNP"="CHR:BP")

write.table(FADS1_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_FADS1_PES_BIP_C_T_cons.txt",
            row.names = F, quote = F)

## FES - Penalised regression (liberal) ## 

FES <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/BIP_TRAIN_lassosum_liberal/BIP_TRAIN_UKBB_LASSO_FES.validate.rds")

FES_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/BIP_TRAIN_lassosum_liberal/BIP_TRAIN_UKBB_LASSO_FES.lassosum.pipeline.rds")

## Summary statistics for FES

FES_weights <- cbind(FES_pipeline$sumstats, FES$best.beta)

FES_NZ <- FES_weights %>% filter(FES$best.beta != 0)

FES_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/Liberal_boundaries_pathway_specific_variants/Liberal_bound_BIP_FES", header = T)

FES_sumstats$order <-  seq.int(nrow(FES_sumstats))

Merged_FES <- merge(FES_weights, FES_sumstats, by = "order")

FES_output <- Merged_FES %>% select("CHR:BP", A1.x, 'FES$best.beta')

FES_output <- rename(FES_output, "A1"="A1.x", "FES_PES"="FES$best.beta", "SNP"="CHR:BP")

write.table(FES_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_FES_PES_BIP_penalised_reg_lib.txt",
            row.names = F, quote = F)

## GRIN2A - Penalised regression (liberal) ## 

GRIN2A <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/BIP_TRAIN_lassosum_liberal/BIP_TRAIN_UKBB_LASSO_GRIN2A.validate.rds")

GRIN2A_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/BIP_TRAIN_lassosum_liberal/BIP_TRAIN_UKBB_LASSO_GRIN2A.lassosum.pipeline.rds")

## Summary statistics for GRIN2A

GRIN2A_weights <- cbind(GRIN2A_pipeline$sumstats, GRIN2A$best.beta)

GRIN2A_NZ <- GRIN2A_weights %>% filter(GRIN2A$best.beta != 0)

GRIN2A_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/Liberal_boundaries_pathway_specific_variants/Liberal_bound_BIP_GRIN2A", header = T)

GRIN2A_sumstats$order <-  seq.int(nrow(GRIN2A_sumstats))

Merged_GRIN2A <- merge(GRIN2A_weights, GRIN2A_sumstats, by = "order")

GRIN2A_output <- Merged_GRIN2A %>% select("CHR:BP", A1.x, 'GRIN2A$best.beta')

GRIN2A_output <- rename(GRIN2A_output, "A1"="A1.x", "GRIN2A_PES"="GRIN2A$best.beta", "SNP"="CHR:BP")

write.table(GRIN2A_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_GRIN2A_PES_BIP_penalised_reg_lib.txt",
            row.names = F, quote = F)

## RPS17 - penalised regression (conservative) ##

RPS17 <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_BIP_UKBB_TRAIN/BIP_TRAIN_UKBB_LASSO_BIP_TRAIN_RPS17_network.validate.rds")

RPS17_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_BIP_UKBB_TRAIN/BIP_TRAIN_UKBB_LASSO_BIP_TRAIN_RPS17_network.lassosum.pipeline.rds")

## Summary statistics for RPS17

RPS17_weights <- cbind(RPS17_pipeline$sumstats, RPS17$best.beta)

RPS17_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/BIP_RPS17_network.txt", header = T)

RPS17_sumstats$order <-  seq.int(nrow(RPS17_sumstats))

Merged_RPS17 <- merge(RPS17_weights, RPS17_sumstats, by = "order")

RPS17_output <- Merged_RPS17 %>% select('CHR:BP', A1.x, 'RPS17$best.beta')

RPS17_output <- rename(RPS17_output, "A1"="A1.x", "RPS17_PES"="RPS17$best.beta", "SNP"="CHR:BP")

write.table(RPS17_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_RPS17_PES_BIP_penalised_reg_cons.txt",
            row.names = F, quote = F)

## PCCB - penalised regression (conservative) ##

PCCB <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_BIP_UKBB_TRAIN/BIP_TRAIN_UKBB_LASSO_BIP_TRAIN_PCCB_network.validate.rds")

PCCB_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_BIP_UKBB_TRAIN/BIP_TRAIN_UKBB_LASSO_BIP_TRAIN_PCCB_network.lassosum.pipeline.rds")

## Summary statistics for PCCB

PCCB_weights <- cbind(PCCB_pipeline$sumstats, PCCB$best.beta)

PCCB_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/BIP_PCCB_network.txt", header = T)

PCCB_sumstats$order <-  seq.int(nrow(PCCB_sumstats))

Merged_PCCB <- merge(PCCB_weights, PCCB_sumstats, by = "order")

PCCB_output <- Merged_PCCB %>% select('CHR:BP', A1.x, 'PCCB$best.beta')

PCCB_output <- rename(PCCB_output, "A1"="A1.x", "PCCB_PES"="PCCB$best.beta", "SNP"="CHR:BP")

write.table(PCCB_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/HCS_cohort_profiling/CHR_BP_HCS_weights/CHR_BP_PCCB_PES_BIP_penalised_reg_cons.txt",
            row.names = F, quote = F)
