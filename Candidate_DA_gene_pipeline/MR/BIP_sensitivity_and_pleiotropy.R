###################################

## Sensitivity and pleiotropy analyses for BIP candidate e-DAG

## William Reay (2020)

###################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/")

library(TwoSampleMR)
library(dplyr)
library(data.table)
library(readxl)
library(coloc)
library(ieugwasr)

## Read in BIP outcome data and format as needed

BIP_RAW <- fread("Raw_summary_statistics/daner_PGC_BIP32b_mds7a_0416a", header = T, fill = T)

BIP_RAW$Beta <- log(BIP_RAW$OR)

## MAP2K2 IV is a trans-effect so will just focus on FADS1 in cortex

Cortex_eQTLs <- fread("MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Final_IVs_LD_clumped_cortex.txt", header = T)

## Disregard id column from LD clumping step

Cortex_eQTLs <- Cortex_eQTLs[, -18]

Cortex_eQTL_IVs <- format_data(dat = Cortex_eQTLs, type = "exposure",
                               effect_allele_col = "A1", other_allele_col = "A2",
                               beta_col = "b", se_col = "SE", phenotype_col = "Probe", snp_col = "rsid", 
                               pval_col = "pval.x")

FADS1_IV <- Cortex_eQTL_IVs %>% filter(exposure == "ENSG00000149485.19")



## Perform colocalisation between genic boundaries and SZ - 1 mb upstream of gene included

BIP_RAW$N <- BIP_RAW$Nca+BIP_RAW$Nco

BIP_RAW$Var <- BIP_RAW$SE^2

## Get all FADS1 SNPs (1 mb flank)

FADS1_all <- fread("MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Coloc_input_FADS1.txt", header = T)

FADS1_all$Var <- FADS1_all$SE^2

Merged_FADS1 <- merge(FADS1_all, BIP_RAW, by="SNP")

Merged_FADS1$samplesize.exposure <- 2907


## Define function to perform colocalisation

Coloc_func <- function(Merged_sumstats, CHR, START, END) {
  
  # Define region in SZ summ stats and filter out SZ info
  BIP_df <- Merged_sumstats %>% filter(CHR == CHR & BP.y > START & BP.y < END) %>% select(SNP, Beta, Var.y, N)
  
  # Define region in eQTL summ stats
  eQTL_df <- Merged_sumstats %>% filter(Chr == CHR & BP.x > START & BP.x < END) %>% select(SNP, b, Var.x, samplesize.exposure)
  
  eQTL_df <- rename(eQTL_df, "N"="samplesize.exposure")
  
  ## Define coloc list input
  
  eQTL=list(beta=eQTL_df$b, varbeta=eQTL_df$Var.x, N=eQTL_df$N, sdY=1, type="quant")
  
  BIP=list(beta=BIP_df$Beta, varbeta=BIP_df$Var.y, N=BIP_df$N,s=0.3922,type="cc")
  
  # Perform colocalisation under default priors - assume SD of 1 for CRP as IRNT
  coloc_df <- coloc.abf(eQTL, BIP, MAF=NULL)
  
  return(coloc_df)
}

FADS1_coloc <- Coloc_func(Merged_FADS1, 11, 60567099, 62584475)

sensitivity(FADS1_coloc, rule="H3 > 0.8")


## FADS1 -  rs174530

rs174530_pheWAS <- phewas("rs174530", pval = 1)

## Flip alleles if effect allele is not A and the non-effect allele is not G

mismatch_FADS1 <- which(rs174530_pheWAS$ea != "A" & rs174530_pheWAS$nea != "G",arr.ind=TRUE)
rs174530_pheWAS[mismatch_FADS1,]$beta <- rs174530_pheWAS[mismatch_FADS1,]$beta*-1
rs174530_pheWAS[mismatch_FADS1,]$ea <- rs174530_pheWAS[mismatch_FADS1,]$nea
rs174530_pheWAS[mismatch_FADS1,]$nea <- rs174530_pheWAS[mismatch_FADS1,]$ea

## Perform phenome-wide pheWAS

rs174530_pheWAS$MR_beta <- rs174530_pheWAS$beta/0.416479
rs174530_pheWAS$MR_SE <- abs(rs174530_pheWAS$se/0.416479)
rs174530_pheWAS$MR_pval <- pnorm(abs(rs174530_pheWAS$MR_beta)/rs174530_pheWAS$MR_SE, lower.tail=FALSE) * 2

write.table(rs174530_pheWAS, file="MR_results/BIP_ALL/BIP_sensitivity_and_pleiotropy/FADS1_SNP_pheWAS_and_MR_pheWAS.txt",
            sep = "\t", row.names = F, quote = F)

## Convert to SMR format

BIP_RAW <- rename(BIP_RAW, "n"="N", "b"="Beta", "se"="SE", "freq"="FRQ_U_31358", "p"="P")

BIP_RAW <- BIP_RAW %>% select(SNP, A1, A2, freq, b, se, p, n)

write.table(BIP_RAW, file="SMR/BIP_SMR.sumstats", sep = "\t", row.names = F, quote = F)

