####################################

## pQTL and eQTL MR - SZ PGC3 (all invidiuals)

## William Reay (2021)

#####################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/")

library(TwoSampleMR)
library(dplyr)
library(data.table)
library(readxl)
library(biomaRt)

## Load MR automated function

load("MR_input/MR_func.Rda")

## Read in SZ PGC3 ALL outcome data and format as needed

SZ_RAW <- fread("Raw_summary_statistics/daner_PGC_SCZ_w3_90_0418b.gz", header = T, fill = T)

SZ_RAW$Beta <- log(SZ_RAW$OR)


###### Blood pQTL - Zheng et al ##########

Blood_pQTL <- read_excel("MR_input/Protein/Blood/Tier_1_Zheng_blood_pQTL_IV_less_than_equal_3_proteins.xlsx")

Blood_pQTL_IVs <- format_data(dat = Blood_pQTL, type = "exposure",
                              effect_allele_col = "Effect_allele", other_allele_col = "Other_allele",
                              beta_col = "Beta", se_col = "SE", pval_col = "P_value", phenotype_col = "Exposure")


## Call function

SZ_blood_pQTL <- Harmonise_and_MR(Blood_pQTL_IVs, SZ_RAW, SZ)

write.table(SZ_blood_pQTL, file="MR_results/SZ_PGC3_ALL/Blood_pQTL_SZ_ALL.txt",
            sep = "\t", row.names = F, quote = F)

###### Brain pQTL - Wingo et al ##########

Raw_brain <- fread("MR_input/Protein/Brain/Processed_raw_input_ROSMAP_consistent_with_Banner_3_or_less_proteins.txt", header = T)

Unfiltered_Brain_pQTL_IVs <- format_data(dat = Raw_brain, type = "exposure",
                              effect_allele_col = "A1", other_allele_col = "A2",
                              beta_col = "Beta", se_col = "SE", pval_col = "P", phenotype_col = "Protein_GeneSymbol")

Filtered_brain_pQTL_IVs <- clump_data(Unfiltered_Brain_pQTL_IVs)

save(Filtered_brain_pQTL_IVs, file="MR_input/Protein/Blood/Wingo_pQTL_LD_clumped.Rda")

SZ_brain_pQTL <- Harmonise_and_MR(Filtered_brain_pQTL_IVs, SZ_RAW, SZ)

write.table(SZ_brain_pQTL, file="MR_results/SZ_PGC3_ALL/Brain_pQTL_SZ_ALL.txt",
            sep = "\t", row.names = F, quote = F)


###### Brain eQTL - Metabrain ##########

## Cortex

Cortex_eQTLs <- fread("MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Final_IVs_LD_clumped_cortex.txt", header = T)

## Disregard id column from LD clumping step

Cortex_eQTLs <- Cortex_eQTLs[, -18]

Gene_names_eQTL <- merge(Cortex_eQTLs, Cortex_eQTLs, by="Gene")

Cortex_eQTL_IVs <- format_data(dat = Cortex_eQTLs, type = "exposure",
                            effect_allele_col = "A1", other_allele_col = "A2",
                            beta_col = "b", se_col = "SE", phenotype_col = "Probe", snp_col = "rsid", 
                            pval_col = "pval.x")

SZ_cortex <- Harmonise_and_MR(Cortex_eQTL_IVs, SZ_RAW, SZ)

write.table(SZ_cortex, file="MR_results/SZ_PGC3_ALL/Cortex_eQTL_SZ_ALL.txt",
            sep = "\t", row.names = F, quote = F)

## Basal ganglia

BG_eQTLs <- fread("MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Final_IVs_LD_clumped_basal_ganglia.txt", header = T)

## Disregard id column from LD clumping step

BG_eQTLs <- BG_eQTLs[, -18]


BG_eQTL_IVs <- format_data(dat = BG_eQTLs, type = "exposure",
                               effect_allele_col = "A1", other_allele_col = "A2",
                               beta_col = "b", se_col = "SE", phenotype_col = "Probe", snp_col = "rsid", 
                               pval_col = "pval.x")

SZ_BG <- Harmonise_and_MR(BG_eQTL_IVs, SZ_RAW, SZ)

write.table(SZ_BG, file="MR_results/SZ_PGC3_ALL/Basal_ganglia_eQTL_SZ_ALL.txt",
            sep = "\t", row.names = F, quote = F)

## Cerebellum

CB_eQTLs <- fread("MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Final_IVs_LD_clumped_cerebellum.txt", header = T)

## Disregard id column from LD clumping step

CB_eQTLs <- CB_eQTLs[, -18]


CB_eQTL_IVs <- format_data(dat = CB_eQTLs, type = "exposure",
                           effect_allele_col = "A1", other_allele_col = "A2",
                           beta_col = "b", se_col = "SE", phenotype_col = "Probe", snp_col = "rsid", 
                           pval_col = "pval.x")

SZ_CB <- Harmonise_and_MR(CB_eQTL_IVs, SZ_RAW, SZ)

write.table(SZ_CB, file="MR_results/SZ_PGC3_ALL/Cerebellum_eQTL_SZ_ALL.txt",
            sep = "\t", row.names = F, quote = F)

## Whole blood - ignore id column from LD prunining

Whole_blood_eQTLs <- fread("MR_input/mRNA/Blood/Final_blood_IVs_LD_clumped.txt", header = T)

Whole_blood_eQTLs <- Whole_blood_eQTLs[, -17]

Blood_eQTL_IVs <- format_data(dat = Whole_blood_eQTLs, type = "exposure",
                           effect_allele_col = "A1", other_allele_col = "A2",
                           beta_col = "b", se_col = "SE", phenotype_col = "Gene", snp_col = "rsid", 
                           pval_col = "pval.x")

SZ_Blood <- Harmonise_and_MR(Blood_eQTL_IVs, SZ_RAW, SZ)

write.table(SZ_Blood, file="MR_results/SZ_PGC3_ALL/Whole_blood_eQTL_SZ_ALL.txt",
            sep = "\t", row.names = F, quote = F)


