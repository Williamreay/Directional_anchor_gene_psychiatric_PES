####################################

## Extracting the score info and output scoring files for best performing PES/PRS

## BIP

## William Reay (2021)

####################################

library(data.table)
library(dplyr)

## PRS - penalised regression ## 

PRS <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_BIP_UKBB_TRAIN/BIP_TRAIN_UKBB_LASSO_BIP_TRAIN_PRS.validate.rds")

PRS_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_BIP_UKBB_TRAIN/BIP_TRAIN_UKBB_LASSO_BIP_TRAIN_PRS.lassosum.pipeline.rds")

## Summary statistics for PRS

PRS_weights <- cbind(PRS_pipeline$sumstats, PRS$best.beta)

BIP_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/QC_sumstats/BIP_QC_sumstats.txt.gz", header = T)

BIP_sumstats$order <-  seq.int(nrow(BIP_sumstats))

Merged_PRS <- merge(PRS_weights, BIP_sumstats, by = "order")

PRS_output <- Merged_PRS %>% select(SNP, A1.x, 'PRS$best.beta')

PRS_output <- rename(PRS_output, "A1"="A1.x", "BIP_PRS"="PRS$best.beta")

write.table(PRS_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/BIP_best_con_or_lib/BIP_SNP_weights_for_validation/Penalised_regression_BIP_gw_PRS.txt",
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

CACNA1C_output <- Merged_CACNA1C %>% select(SNP, A1.x, 'CACNA1C$best.beta')

CACNA1C_output <- rename(CACNA1C_output, "A1"="A1.x", "CACNA1C_PES"="CACNA1C$best.beta")

write.table(CACNA1C_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/BIP_best_con_or_lib/BIP_SNP_weights_for_validation/CACNA1C_PES_BIP_penalised_reg_lib.txt",
            row.names = F, quote = F)

## FADS1 - C+T (conservative)

FADS1_SNPs_to_include <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/BIP_TRAIN_UKBB/BIP_TRAIN_UKBB_BIP_TRAIN_FADS1_network.snp", header = T)

FADS1_sum_stats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/BIP_FADS1_network.txt", header = T)

FADS1_sum_stats$Beta <- log(FADS1_sum_stats$OR)

Merged_FADS1 <- merge(FADS1_sum_stats, FADS1_SNPs_to_include, by = "SNP")

FADS1_output <- Merged_FADS1 %>% select(SNP, A1, Beta)

FADS1_output <- rename(FADS1_output, "FADS1_PES"="Beta")

write.table(FADS1_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/BIP_best_con_or_lib/BIP_SNP_weights_for_validation/FADS1_PES_BIP_C_T_cons.txt",
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

FES_output <- Merged_FES %>% select(SNP, A1.x, 'FES$best.beta')

FES_output <- rename(FES_output, "A1"="A1.x", "FES_PES"="FES$best.beta")

write.table(FES_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/BIP_best_con_or_lib/BIP_SNP_weights_for_validation/FES_PES_BIP_penalised_reg_lib.txt",
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

GRIN2A_output <- Merged_GRIN2A %>% select(SNP, A1.x, 'GRIN2A$best.beta')

GRIN2A_output <- rename(GRIN2A_output, "A1"="A1.x", "GRIN2A_PES"="GRIN2A$best.beta")

write.table(GRIN2A_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/BIP_best_con_or_lib/BIP_SNP_weights_for_validation/GRIN2A_PES_BIP_penalised_reg_lib.txt",
            row.names = F, quote = F)

## RPS17 - penalised regression (liberal) ##

RPS17 <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/BIP_TRAIN_lassosum_liberal/BIP_TRAIN_UKBB_LASSO_RPS17.validate.rds")

RPS17_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/BIP_TRAIN_lassosum_liberal/BIP_TRAIN_UKBB_LASSO_RPS17.lassosum.pipeline.rds")

## Summary statistics for RPS17

RPS17_weights <- cbind(RPS17_pipeline$sumstats, RPS17$best.beta)

RPS17_NZ <- RPS17_weights %>% filter(RPS17$best.beta != 0)

RPS17_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/Liberal_boundaries_pathway_specific_variants/Liberal_bound_BIP_RPS17", header = T)

RPS17_sumstats$order <-  seq.int(nrow(RPS17_sumstats))

Merged_RPS17 <- merge(RPS17_weights, RPS17_sumstats, by = "order")

RPS17_output <- Merged_RPS17 %>% select(SNP, A1.x, 'RPS17$best.beta')

RPS17_output <- rename(RPS17_output, "A1"="A1.x", "RPS17_PES"="RPS17$best.beta")

write.table(RPS17_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/BIP_best_con_or_lib/BIP_SNP_weights_for_validation/RPS17_PES_BIP_penalised_reg_cons.txt",
            row.names = F, quote = F)

## PCCB - penalised regression (liberal) ##

PCCB <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/BIP_TRAIN_lassosum_liberal/BIP_TRAIN_UKBB_LASSO_PCCB.validate.rds")

PCCB_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/BIP_TRAIN_lassosum_liberal/BIP_TRAIN_UKBB_LASSO_PCCB.lassosum.pipeline.rds")

## Summary statistics for PCCB

PCCB_weights <- cbind(PCCB_pipeline$sumstats, PCCB$best.beta)

PCCB_NZ <- PCCB_weights %>% filter(`PCCB$best.beta` != 0)

PCCB_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/Liberal_boundaries_pathway_specific_variants/Liberal_bound_BIP_PCCB", header = T)

PCCB_sumstats$order <-  seq.int(nrow(PCCB_sumstats))

Merged_PCCB <- merge(PCCB_weights, PCCB_sumstats, by = "order")

PCCB_output <- Merged_PCCB %>% select(SNP, A1.x, 'PCCB$best.beta')

PCCB_output <- rename(PCCB_output, "A1"="A1.x", "PCCB_PES"="PCCB$best.beta")

write.table(PCCB_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/BIP_best_con_or_lib/BIP_SNP_weights_for_validation/PCCB_PES_BIP_penalised_reg_cons.txt",
            row.names = F, quote = F)

PCBB_NZ <- PCCB_weights %>% filter(PCCB$best.beta != 0)
