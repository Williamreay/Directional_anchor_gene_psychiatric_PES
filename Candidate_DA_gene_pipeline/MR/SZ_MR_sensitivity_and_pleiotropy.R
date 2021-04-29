###################################

## Sensitivity and pleiotropy analyses for SZ candidate e-DAG

## William Reay (2020)

###################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/")

library(TwoSampleMR)
library(dplyr)
library(data.table)
library(readxl)
library(coloc)
library(ieugwasr)

## Read in SZ PGC3 ALL outcome data and format as needed

SZ_RAW <- fread("Raw_summary_statistics/daner_PGC_SCZ_w3_90_0418b.gz", header = T, fill = T)

SZ_RAW$Beta <- log(SZ_RAW$OR)

## Genes to consider:
## Blood - NEK1 and PTK2B
## Cortex - PCCB

## Blood ##

Whole_blood_eQTLs <- fread("MR_input/mRNA/Blood/Final_blood_IVs_LD_clumped.txt", header = T)

Whole_blood_eQTLs <- Whole_blood_eQTLs[, -17]

Blood_eQTL_IVs <- format_data(dat = Whole_blood_eQTLs, type = "exposure",
                              effect_allele_col = "A1", other_allele_col = "A2",
                              beta_col = "b", se_col = "SE", phenotype_col = "Probe", snp_col = "rsid", 
                              pval_col = "pval.x")

## Blood - select out IVs for genes in question

Genes_to_test <- Blood_eQTL_IVs %>% filter(exposure == "ENSG00000120899" |
                                             exposure == "ENSG00000137601")



## Extract from SZ and harmonise

SZ_Extract <- format_data(dat = SZ_RAW, type = "outcome",
                       snps = Genes_to_test$SNP, snp_col = "SNP", beta_col = "Beta",
                       se_col = "SE", effect_allele_col = "A1",
                       other_allele_col = "A2", pval_col = "P")

Sensitivity_analyses_harmonised <- harmonise_data(Genes_to_test, SZ_Extract, action = 3)

Sensitivity_analyses_harmonised <- Sensitivity_analyses_harmonised %>% filter(mr_keep == TRUE)

write.table(Sensitivity_analyses_harmonised, file="MR_results/SZ_PGC3_ALL/MR_sensitivity_and_pleiotropy_analyses/Blood_Harmonised_df_sensitivity_analyses.txt",
            sep = "\t", row.names = F, quote = F)

## Test reverse causality of SZ on mRNA

Sensitivity_analyses_harmonised$Reverse_beta <- Sensitivity_analyses_harmonised$beta.exposure/Sensitivity_analyses_harmonised$beta.outcome

Sensitivity_analyses_harmonised$Reverse_SE <- abs(Sensitivity_analyses_harmonised$se.exposure/Sensitivity_analyses_harmonised$beta.outcome)

Sensitivity_analyses_harmonised$Reverse_P <- pnorm(abs(Sensitivity_analyses_harmonised$Reverse_beta)/Sensitivity_analyses_harmonised$Reverse_SE, lower.tail=FALSE) * 2

## Perform colocalisation between genic boundaries and SZ - 1 mb upstream of gene included

SZ_RAW$N <- SZ_RAW$Nca+SZ_RAW$Nco

SZ_RAW$Var <- SZ_RAW$SE^2

## Get all SNPs in region from IEUGWAS

## PTK2B

PTK2B_all <- extract_instruments("eqtl-a-ENSG00000120899", p1 = 1, clump = F)

PTK2B_all$Var <- PTK2B_all$se.exposure^2

Merged_PT2KB <- merge(PTK2B_all, SZ_RAW, by="SNP")


## Define function to perform colocalisation

Coloc_func <- function(Merged_sumstats, CHR, START, END) {
  
  # Define region in SZ summ stats and filter out SZ info
  SZ_df <- Merged_sumstats %>% filter(CHR == CHR & BP > START & BP < END) %>% select(SNP, Beta, Var.y, N)
  
  # Define region in eQTL summ stats
  eQTL_df <- Merged_sumstats %>% filter(chr.exposure == CHR & pos.exposure > START & pos.exposure < END) %>% select(SNP, beta.exposure, Var.x, samplesize.exposure)
  
  eQTL_df <- rename(eQTL_df, "N"="samplesize.exposure")
  
  ## Define coloc list input
  
  eQTL=list(beta=eQTL_df$beta.exposure, varbeta=eQTL_df$Var.x, N=eQTL_df$N, sdY=1, type="quant")
  
  SZ=list(beta=SZ_df$Beta, varbeta=SZ_df$Var.y, N=SZ_df$N,s=0.409,type="cc")
  
  # Perform colocalisation under default priors - assume SD of 1 for CRP as IRNT
  coloc_df <- coloc.abf(eQTL, SZ, MAF=NULL)
  
  return(coloc_df)
}

PT2KB_coloc <- Coloc_func(Merged_PT2KB, 8, 26168995,  28316908)

sensitivity(PT2KB_coloc, rule="H3 > 0.8")

## NEK1

NEK1_all <- extract_instruments("eqtl-a-ENSG00000137601", p1 = 1, clump = F)

NEK1_all$Var <- NEK1_all$se.exposure^2

Merged_NEK1 <- merge(NEK1_all, SZ_RAW, by="SNP")

NEK1_coloc <- Coloc_func(Merged_NEK1, 4, 169313960, 171533734)

sensitivity(NEK1_coloc, rule="H3 > 0.8")

## Brain

## PCCB

## Read in raw SNPs on chromosome 3 - average number of samples for cis-eQTLs in PCCB region = 2507

Chromsome_3_MetaBrain_cortex <- fread("MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/PCCB_cortex_chrosome_3.txt", header = T)

Chromsome_3_MetaBrain_cortex$Var <- Chromsome_3_MetaBrain_cortex$SE^2

Chromsome_3_MetaBrain_cortex <- rename(Chromsome_3_MetaBrain_cortex,
                                       "beta.exposure"="b", "chr.exposure"="Chr",
                                       "pos.exposure"="BP")

Chromsome_3_MetaBrain_cortex$samplesize.exposure <- 2970

Merged_PCCB <- merge(Chromsome_3_MetaBrain_cortex, SZ_RAW, by="SNP")

PCCB_coloc <- Coloc_func(Merged_PCCB, 3, 134969182, 137049011)

sensitivity(PCCB_coloc, rule="H3 > 0.8")

## Pleiotropy based pheWAS using IEUGWASR

## SERPINC1 - rs9425780

rs9425780_pheWAS <- phewas("rs9425780", pval = 1)

## Flip alleles if effect allele is not A and the non-effect allele is not G

mismatch_SERPINC1 <- which(rs9425780_pheWAS$ea != "A" & rs9425780_pheWAS$nea != "G",arr.ind=TRUE)
rs9425780_pheWAS[mismatch_SERPINC1,]$beta <- rs9425780_pheWAS[mismatch_SERPINC1,]$beta*-1
rs9425780_pheWAS[mismatch_SERPINC1,]$ea <- rs9425780_pheWAS[mismatch_SERPINC1,]$nea
rs9425780_pheWAS[mismatch_SERPINC1,]$nea <- rs9425780_pheWAS[mismatch_SERPINC1,]$ea

## Perform phenome-wide pheWAS

rs9425780_pheWAS$MR_beta <- rs9425780_pheWAS$beta/0.600407
rs9425780_pheWAS$MR_SE <- abs(rs9425780_pheWAS$se/0.600407)
rs9425780_pheWAS$MR_pval <- pnorm(abs(rs9425780_pheWAS$MR_beta)/rs9425780_pheWAS$MR_SE, lower.tail=FALSE) * 2

write.table(rs9425780_pheWAS, file="MR_results/SZ_PGC3_ALL/MR_sensitivity_and_pleiotropy_analyses/pheWAS/SERPINC1_SNP_pheWAS_and_MR_pheWAS.txt",
            sep = "\t", row.names = F, quote = F)

## MAPK3 - rs4548895 (main SNP from LOO)

rs4548895_pheWAS <- phewas("rs4548895", pval = 1)

## Flip alleles if effect allele is not A and the non-effect allele is not G

mismatch_MAPK3 <- which(rs4548895_pheWAS$ea != "A" & rs4548895_pheWAS$nea != "G",arr.ind=TRUE)

## No mismatches

## Perform phenome-wide MR

rs4548895_pheWAS$MR_beta <- rs4548895_pheWAS$beta/-0.560603
rs4548895_pheWAS$MR_SE <- abs(rs4548895_pheWAS$se/-0.560603)
rs4548895_pheWAS$MR_pval <- pnorm(abs(rs4548895_pheWAS$MR_beta)/rs4548895_pheWAS$MR_SE, lower.tail=FALSE) * 2

write.table(rs4548895_pheWAS, file="MR_results/SZ_PGC3_ALL/MR_sensitivity_and_pleiotropy_analyses/pheWAS/MAPK3_rs4548895_SNP_pheWAS_and_MR_pheWAS.txt",
            sep = "\t", row.names = F, quote = F)

## NEK1 - rs7696352 as the SNP with larger effect

rs55909907_pheWAS <- phewas("rs55909907", pval = 1)

## Flip alleles if effect allele is not A and the non-effect allele is not G

mismatch_NEK1 <- which(rs55909907_pheWAS$ea != "T" & rs55909907_pheWAS$nea != "C",arr.ind=TRUE)

rs55909907_pheWAS[mismatch_NEK1,]$beta <- rs55909907_pheWAS[mismatch_NEK1,]$beta*-1
rs55909907_pheWAS[mismatch_NEK1,]$ea <- rs55909907_pheWAS[mismatch_NEK1,]$nea
rs55909907_pheWAS[mismatch_NEK1,]$nea <- rs55909907_pheWAS[mismatch_NEK1,]$ea

## Perform phenome-wide MR

rs55909907_pheWAS$MR_beta <- rs55909907_pheWAS$beta/-0.147534
rs55909907_pheWAS$MR_SE <- abs(rs55909907_pheWAS$se/-0.147534)
rs55909907_pheWAS$MR_pval <- pnorm(abs(rs55909907_pheWAS$MR_beta)/rs55909907_pheWAS$MR_SE, lower.tail=FALSE) * 2

write.table(rs55909907_pheWAS, file="MR_results/SZ_PGC3_ALL/MR_sensitivity_and_pleiotropy_analyses/pheWAS/NEK1_rs55909907_SNP_pheWAS_and_MR_pheWAS.txt",
            sep = "\t", row.names = F, quote = F)

## PT2K2B

rs3905272_pheWAS <- phewas("rs3905272", pval = 1)

## Flip alleles if effect allele is not A and the non-effect allele is not G

mismatch_PTK2B <- which(rs3905272_pheWAS$ea != "C" & rs3905272_pheWAS$nea != "T",arr.ind=TRUE)

rs3905272_pheWAS[mismatch_PTK2B,]$beta <- rs3905272_pheWAS[mismatch_PTK2B,]$beta*-1
rs3905272_pheWAS[mismatch_PTK2B,]$ea <- rs3905272_pheWAS[mismatch_PTK2B,]$nea
rs3905272_pheWAS[mismatch_PTK2B,]$nea <- rs3905272_pheWAS[mismatch_PTK2B,]$ea

## Perform phenome-wide MR

rs3905272_pheWAS$MR_beta <- rs3905272_pheWAS$beta/0.236684
rs3905272_pheWAS$MR_SE <- abs(rs3905272_pheWAS$se/0.236684)
rs3905272_pheWAS$MR_pval <- pnorm(abs(rs3905272_pheWAS$MR_beta)/rs3905272_pheWAS$MR_SE, lower.tail=FALSE) * 2

write.table(rs3905272_pheWAS, file="MR_results/SZ_PGC3_ALL/MR_sensitivity_and_pleiotropy_analyses/pheWAS/PTK2B_rs3905272_SNP_pheWAS_and_MR_pheWAS.txt",
            sep = "\t", row.names = F, quote = F)

## PCCB (cortex effect size pheWAS)

## Cortex

Cortex_eQTLs <- fread("MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Final_IVs_LD_clumped_cortex.txt", header = T)

## Disregard id column from LD clumping step

Cortex_eQTLs <- Cortex_eQTLs[, -18]

Cortex_eQTL_IVs <- format_data(dat = Cortex_eQTLs, type = "exposure",
                               effect_allele_col = "A1", other_allele_col = "A2",
                               beta_col = "b", se_col = "SE", phenotype_col = "Probe", snp_col = "rsid", 
                               pval_col = "pval.x")

PCCB_IV <- Cortex_eQTL_IVs %>% filter(exposure == "ENSG00000114054.14")

## PCCB IV pheWAS - rs480330  

rs480330_pheWAS <- phewas("rs480330", pval = 1)

## Flip alleles if effect allele is not A and the non-effect allele is not G

mismatch_PCCB <- which(rs480330_pheWAS$ea != "C" &  rs480330_pheWAS$nea != "T",arr.ind=TRUE)

rs480330_pheWAS[mismatch_PCCB,]$beta <- rs480330_pheWAS[mismatch_PCCB,]$beta*-1
rs480330_pheWAS[mismatch_PCCB,]$ea <- rs480330_pheWAS[mismatch_PCCB,]$nea
rs480330_pheWAS[mismatch_PCCB,]$nea <- rs480330_pheWAS[mismatch_PCCB,]$ea

## Perform phenome-wide MR

rs480330_pheWAS$MR_beta <- rs480330_pheWAS$beta/-0.325655
rs480330_pheWAS$MR_SE <- abs(rs480330_pheWAS$se/-0.325655)
rs480330_pheWAS$MR_pval <- pnorm(abs(rs480330_pheWAS$MR_beta)/rs480330_pheWAS$MR_SE, lower.tail=FALSE) * 2

write.table(rs480330_pheWAS, file="MR_results/SZ_PGC3_ALL/MR_sensitivity_and_pleiotropy_analyses/pheWAS/PCCB_rs480330_SNP_pheWAS_and_MR_pheWAS.txt",
            sep = "\t", row.names = F, quote = F)

## Format to SMR format

SZ_RAW <- rename(SZ_RAW, "n"="N", "b"="Beta", "se"="SE", "freq"="FRQ_U_94015")

SZ_RAW <- rename(SZ_RAW, "p"="P")

SZ_RAW <- SZ_RAW %>% select(SNP, A1, A2, freq, b, se, p, n)

write.table(SZ_RAW, file="SMR/SZ_PGC3_SMR.sumstats", sep = "\t", row.names = F, quote = F)
