####################################

## pQTL and eQTL MR - Bipolar (all invidiuals)

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

BIP_RAW <- fread("Raw_summary_statistics/daner_PGC_BIP32b_mds7a_0416a", header = T, fill = T)

BIP_RAW$Beta <- log(BIP_RAW$OR)


###### Blood pQTL - Zheng et al ##########

Blood_pQTL <- read_excel("MR_input/Protein/Blood/Tier_1_Zheng_blood_pQTL_IV_less_than_equal_3_proteins.xlsx")

Blood_pQTL_IVs <- format_data(dat = Blood_pQTL, type = "exposure",
                              effect_allele_col = "Effect_allele", other_allele_col = "Other_allele",
                              beta_col = "Beta", se_col = "SE", pval_col = "P_value", phenotype_col = "Exposure")

## Remove rs67136035 SNP with undefined other allele

Blood_pQTL_IVs <- Blood_pQTL_IVs %>% filter(SNP != "rs67136035")

## Call function

BIP_blood_pQTL <- Harmonise_and_MR(Blood_pQTL_IVs, BIP_RAW, BIP)

write.table(BIP_blood_pQTL, file="MR_results/BIP_ALL/Blood_pQTL_BIP_ALL.txt",
            sep = "\t", row.names = F, quote = F)

###### Brain pQTL - Wingo et al ##########

## Open .Rdat file

load("MR_input/Protein/Brain/Wingo_pQTL_LD_clumped.Rda")

BIP_brain_pQTL <- Harmonise_and_MR(Filtered_brain_pQTL_IVs, BIP_RAW, BIP)

write.table(BIP_brain_pQTL, file="MR_results/BIP_ALL/Brain_pQTL_BIP_ALL.txt",
            sep = "\t", row.names = F, quote = F)

## FDR correction

pQTL_BIP <- rbind(BIP_blood_pQTL, BIP_brain_pQTL)
pQTL_BIP$FWER <- p.adjust(pQTL_BIP$pval, method = "bonferroni")
pQTL_BIP$FDR <- p.adjust(pQTL_BIP$pval, method = "fdr")

write.table(pQTL_BIP, file="MR_results/BIP_ALL/pQTL_blood_brain_BIP_FWER_FDR.txt",
            sep = "\t", row.names = F, quote = F)


###### Brain eQTL - Metabrain ##########

## Cortex

Cortex_eQTLs <- fread("MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Final_IVs_LD_clumped_cortex.txt", header = T)

## Disregard id column from LD clumping step

Cortex_eQTLs <- Cortex_eQTLs[, -18]

Cortex_eQTL_IVs <- format_data(dat = Cortex_eQTLs, type = "exposure",
                            effect_allele_col = "A1", other_allele_col = "A2",
                            beta_col = "b", se_col = "SE", phenotype_col = "Probe", snp_col = "rsid", 
                            pval_col = "pval.x")

BIP_cortex <- Harmonise_and_MR(Cortex_eQTL_IVs, BIP_RAW, BIP)

write.table(BIP_cortex, file="MR_results/BIP_ALL/Cortex_eQTL_BIP_ALL.txt",
            sep = "\t", row.names = F, quote = F)

## Basal ganglia

BG_eQTLs <- fread("MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Final_IVs_LD_clumped_basal_ganglia.txt", header = T)

## Disregard id column from LD clumping step

BG_eQTLs <- BG_eQTLs[, -18]


BG_eQTL_IVs <- format_data(dat = BG_eQTLs, type = "exposure",
                               effect_allele_col = "A1", other_allele_col = "A2",
                               beta_col = "b", se_col = "SE", phenotype_col = "Probe", snp_col = "rsid", 
                               pval_col = "pval.x")

BIP_BG <- Harmonise_and_MR(BG_eQTL_IVs, BIP_RAW, BIP)

write.table(BIP_BG, file="MR_results/BIP_ALL/Basal_ganglia_eQTL_BIP_ALL.txt",
            sep = "\t", row.names = F, quote = F)

## Cerebellum

CB_eQTLs <- fread("MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Final_IVs_LD_clumped_cerebellum.txt", header = T)

## Disregard id column from LD clumping step

CB_eQTLs <- CB_eQTLs[, -18]


CB_eQTL_IVs <- format_data(dat = CB_eQTLs, type = "exposure",
                           effect_allele_col = "A1", other_allele_col = "A2",
                           beta_col = "b", se_col = "SE", phenotype_col = "Probe", snp_col = "rsid", 
                           pval_col = "pval.x")

BIP_CB <- Harmonise_and_MR(CB_eQTL_IVs, BIP_RAW, BIP)

write.table(BIP_CB, file="MR_results/BIP_ALL/Cerebellum_eQTL_BIP_ALL.txt",
            sep = "\t", row.names = F, quote = F)

## Whole blood - ignore id column from LD prunining

Whole_blood_eQTLs <- fread("MR_input/mRNA/Blood/Final_blood_IVs_LD_clumped.txt", header = T)

Whole_blood_eQTLs <- Whole_blood_eQTLs[, -17]

Blood_eQTL_IVs <- format_data(dat = Whole_blood_eQTLs, type = "exposure",
                           effect_allele_col = "A1", other_allele_col = "A2",
                           beta_col = "b", se_col = "SE", phenotype_col = "Gene", snp_col = "rsid", 
                           pval_col = "pval.x")

BIP_Blood <- Harmonise_and_MR(Blood_eQTL_IVs, BIP_RAW, BIP)

write.table(BIP_Blood, file="MR_results/BIP_ALL/Whole_blood_eQTL_BIP_ALL.txt",
            sep = "\t", row.names = F, quote = F)


