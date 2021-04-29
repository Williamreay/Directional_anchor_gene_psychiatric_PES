########################################

## LD clumping the MetaBrain eQTL SNPs to select IVs

## William Reay (2021)

#########################################

library(dplyr)
library(data.table)
library(ieugwasr)
library(optparse)

MetaBrain_raw <- fread("~/data/users/william/MetaBrain_eQTL_SMR_format/MetaBrain_Cortex_less_than_4_genes_eQTL_GW_sig.txt", header = T)

MetaBrain_raw <- rename(MetaBrain_raw, "rsid"="SNP", "pval"="p")

List_genes <- as.list(MetaBrain_raw$Gene)

LD_clump <- function(Exposure, Gene_name) {
  eQTL <- Exposure %>% filter(Gene == Gene_name)
  Gene <- ld_clump(dplyr::tibble(rsid=eQTL$rsid, pval=eQTL$pval), 
                   clump_r2 =  0.001, clump_p = 5e-08, plink_bin = "/home/control/R/x86_64-pc-linux-gnu-library/4.0/genetics.binaRies/bin/plink_Linux",
                   bfile = "/home/control/data/genomic_data/1000g_LDREF/EUR")
  return(Gene)
}


Cortex <- ld_clump(dplyr::tibble(rsid=MetaBrain_raw$rsid, pval=MetaBrain_raw$pval), 
                   clump_r2 =  0.001, clump_p = 5e-08, plink_bin = "/home/control/R/x86_64-pc-linux-gnu-library/4.0/genetics.binaRies/bin/plink_Linux",
                   bfile = "/home/control/data/genomic_data/1000g_LDREF/EUR")

Merge_cortex <- merge(MetaBrain_raw, Cortex, by="rsid")

## Remove non-matching p value SNPs, i.e. duplicated for that exposure

Final_merged_cortex_IVs <- Merge_cortex %>% filter(pval.x == pval.y)

write.table(Final_merged_cortex_IVs, file="data/users/william/MetaBrain_eQTL_SMR_format/Final_IVs_LD_clumped_cortex.txt",
            sep = "\t", row.names = F, quote = F)



## Basal Ganglia 

BG_MetaBrain_raw <- fread("~/data/users/william/MetaBrain_eQTL_SMR_format/Basal_ganglia_Less_than_4_genes_per_QTL_MetaBrain.txt", header = T)

BG_MetaBrain_raw <- rename(BG_MetaBrain_raw, "pval"="p")

BG <- ld_clump(dplyr::tibble(rsid=BG_MetaBrain_raw$rsid, pval=BG_MetaBrain_raw$pval), 
               clump_r2 =  0.001, clump_p = 5e-08, plink_bin = "/home/control/R/x86_64-pc-linux-gnu-library/4.0/genetics.binaRies/bin/plink_Linux",
               bfile = "/home/control/data/genomic_data/1000g_LDREF/EUR")

Merge_BG <- merge(BG_MetaBrain_raw, BG, by="rsid")

## Remove non-matching p value SNPs, i.e. duplicated for that exposure

Final_merged_BG_IVs <- Merge_BG %>% filter(pval.x == pval.y)

write.table(Final_merged_BG_IVs, file="data/users/william/MetaBrain_eQTL_SMR_format/Final_IVs_LD_clumped_basal_ganglia.txt",
            sep = "\t", row.names = F, quote = F)

## Cerebellum

CB_MetaBrain_raw <- fread("~/data/users/william/MetaBrain_eQTL_SMR_format/Cerebellum_Less_than_4_genes_per_QTL_MetaBrain.txt", header = T)

CB_MetaBrain_raw <- rename(CB_MetaBrain_raw, "pval"="p", "rsid"="SNP")

CB <- ld_clump(dplyr::tibble(rsid=CB_MetaBrain_raw$rsid, pval=CB_MetaBrain_raw$pval), 
               clump_r2 =  0.001, clump_p = 5e-08, plink_bin = "/home/control/R/x86_64-pc-linux-gnu-library/4.0/genetics.binaRies/bin/plink_Linux",
               bfile = "/home/control/data/genomic_data/1000g_LDREF/EUR")

Merge_CB <- merge(CB_MetaBrain_raw, CB, by="rsid")

## Remove non-matching p value SNPs, i.e. duplicated for that exposure

Final_merged_CB_IVs <- Merge_CB %>% filter(pval.x == pval.y)

write.table(Final_merged_CB_IVs, file="data/users/william/MetaBrain_eQTL_SMR_format/Final_IVs_LD_clumped_cerebellum.txt",
            sep = "\t", row.names = F, quote = F)

## Blood eQTL gen

Blood <- fread("~/data/users/william/MetaBrain_eQTL_SMR_format/RAW_no_LD_pruning_Blood_IVs_three_or_less_genes.txt", header = T)

Blood <- rename(Blood, "pval"="p", "rsid"="SNP")

Blood_clump <- ld_clump(dplyr::tibble(rsid=Blood$rsid, pval=Blood$pval), 
                        clump_r2 =  0.001, clump_p = 5e-08, plink_bin = "/home/control/R/x86_64-pc-linux-gnu-library/4.0/genetics.binaRies/bin/plink_Linux",
                        bfile = "/home/control/data/genomic_data/1000g_LDREF/EUR")

Merge_blood <- merge(Blood, Blood_clump, by="rsid")

## Remove non-matching p value SNPs, i.e. duplicated for that exposure

Final_merged_blood_IVs <- Merge_blood %>% filter(pval.x == pval.y)


write.table(Final_merged_blood_IVs, file="data/users/william/MetaBrain_eQTL_SMR_format/Final_blood_IVs_LD_clumped.txt",
            sep = "\t", row.names = F, quote = F)