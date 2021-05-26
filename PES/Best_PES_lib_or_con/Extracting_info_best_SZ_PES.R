####################################

## Extracting the score info and output scoring files for best performing PES/PRS

## SZ

## William Reay (2021)

####################################

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

PRS_output <- Merged_PRS %>% select(SNP, A1.x, 'PRS$best.beta')

PRS_output <- rename(PRS_output, "A1"="A1.x", "SZ_PRS"="PRS$best.beta")

write.table(PRS_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/SZ_best_con_or_lib/SZ_SNP_weights_for_validation/SZ_PRS_penalised_reg_weights.txt",
            row.names = F, quote = F)

## CACNA1C - C+T (conservative boundaries) ## 

CACNA1C_SNPs_to_include <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/UKBB_SZ/SZ_UKBB_SZ_PGC3_ALL_CACNA1C_network.snp", header = T)

CACNA1C_SNPs_to_include <- CACNA1C_SNPs_to_include %>% filter(P < 0.005)

CACNA1C_sum_stats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/SZ_PGC3_ALL_CACNA1C_network.txt", header = T)

CACNA1C_sum_stats$Beta <- log(CACNA1C_sum_stats$OR)

Merged_CACNA1C <- merge(CACNA1C_sum_stats, CACNA1C_SNPs_to_include, by = "SNP")

CACNA1C_output <- Merged_CACNA1C %>% select(SNP, A1, Beta)

CACNA1C_output <- rename(CACNA1C_output, "CACNA1C_PES"="Beta")

write.table(CACNA1C_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/SZ_best_con_or_lib/SZ_SNP_weights_for_validation/CACNA1C_PES_SZ_cons_C_T.txt",
            row.names = F, quote = F)

## FADS1 - C+T (conservative boundaries) ## 

FADS1_SNPs_to_include <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/UKBB_SZ/SZ_UKBB_SZ_PGC3_ALL_FADS1_network.snp", header = T)

FADS1_SNPs_to_include <- FADS1_SNPs_to_include %>% filter(P < 0.05)

FADS1_sum_stats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/SZ_PGC3_ALL_FADS1_network.txt", header = T)

FADS1_sum_stats$Beta <- log(FADS1_sum_stats$OR)

Merged_FADS1 <- merge(FADS1_sum_stats, FADS1_SNPs_to_include, by = "SNP")

FADS1_output <- Merged_FADS1 %>% select(SNP, A1, Beta)

FADS1_output <- rename(FADS1_output, "FADS1_PES"="Beta")

write.table(FADS1_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/SZ_best_con_or_lib/SZ_SNP_weights_for_validation/FADS1_PES_SZ_cons_C_T.txt",
            row.names = F, quote = F)

## FES - C+T (liberal) ##

FES_SNPs_to_include <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/SZ_UKBB_liberal/SZ_UKBB_FES.snp", header = T)

FES_SNPs_to_include <- FES_SNPs_to_include %>% filter(P < 0.05)

FES_sum_stats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Liberal_genic_boundaries_DA_gene_PES/Liberal_boundaries_pathway_specific_variants/Liberal_bound_SZ_all_PGC3_FES", header = T)

FES_sum_stats$Beta <- log(FES_sum_stats$OR)

Merged_FES <- merge(FES_sum_stats, FES_SNPs_to_include, by = "SNP")

FES_output <- Merged_FES %>% select(SNP, A1, Beta)

FES_output <- rename(FES_output, "FES_PES"="Beta")

write.table(FES_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/SZ_best_con_or_lib/SZ_SNP_weights_for_validation/FES_PES_SZ_cons_C_T.txt",
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

GRIN2A_output <- Merged_GRIN2A %>% select(SNP, A1.x, 'GRIN2A$best.beta')

GRIN2A_output <- rename(GRIN2A_output, "GRIN2A_PES"="GRIN2A$best.beta", "A1"="A1.x")

write.table(GRIN2A_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/SZ_best_con_or_lib/SZ_SNP_weights_for_validation/GRIN2A_PES_penalised_reg_weights.txt",
            row.names = F, quote = F)

## PCCB - penalised regression ##

PCCB <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_PCCB_network.validate.rds")

PCCB_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_PCCB_network.lassosum.pipeline.rds")

## Summary statistics for PRS

PCCB_weights <- cbind(PCCB_pipeline$sumstats, PCCB$best.beta)

PCCB_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/SZ_PGC3_ALL_PCCB_network.txt", header = T)

PCCB_sumstats$order <-  seq.int(nrow(PCCB_sumstats))

Merged_PCCB <- merge(PCCB_weights, PCCB_sumstats, by = "order")

PCCB_output <- Merged_PCCB %>% select(SNP, A1.x, 'PCCB$best.beta')

PCCB_output <- rename(PCCB_output, "PCCB_PES"="PCCB$best.beta", "A1"="A1.x")

write.table(PCCB_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/SZ_best_con_or_lib/SZ_SNP_weights_for_validation/PCCB_PES_penalised_reg_weights.txt",
            row.names = F, quote = F)

## RPS17 - penalised regression (conservative) ##

RPS17 <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_RPS17_network.validate.rds")

RPS17_pipeline <- readRDS("~/Desktop/SZ_PES_mtCOJO_norm/PES/Lassosum_SZ_UKBB/SZ_UKBB_LASSO_SZ_PGC3_ALL_RPS17_network.lassosum.pipeline.rds")

## Summary statistics for PRS

RPS17_weights <- cbind(RPS17_pipeline$sumstats, RPS17$best.beta)

RPS17_sumstats <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/Pathway_specific_variants/SZ_PGC3_ALL_RPS17_network.txt", header = T)

RPS17_sumstats$order <-  seq.int(nrow(RPS17_sumstats))

Merged_RPS17 <- merge(RPS17_weights, RPS17_sumstats, by = "order")

RPS17_output <- Merged_RPS17 %>% select(SNP, A1.x, 'RPS17$best.beta')

RPS17_output <- rename(RPS17_output, "RPS17_PES"="RPS17$best.beta", "A1"="A1.x")

write.table(RPS17_output, file="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/SZ_best_con_or_lib/SZ_SNP_weights_for_validation/RPS17_PES_penalised_reg_weights.txt",
            row.names = F, quote = F)

