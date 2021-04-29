##########################################

## Processing QTLs for MR input

## William Reay (2021)

##########################################

library(data.table)
library(dplyr)


## Brain pQTLs
## Identify i) GW. sig SNPs in ROSMAP which are also P < 0.05 in the Banner study, ii) Associated with <= 3 proteins

ROSMAP <- fread("~/Desktop/SZ_PES_mtCOJO_norm/MR_input/Protein/Brain/Uniq_ID_ROSMAP.txt", header = T)

GW_sig_ROSMAP <- ROSMAP %>% filter(P < 5e-08)

Banner <- fread("~/Desktop/SZ_PES_mtCOJO_norm/MR_input/Protein/Brain/Formatted_Banner_pQTL.tab.txt", header = T)

Nom_sig_Banner <- Banner %>% filter(P < 0.05)

Nom_sig_Banner <- Nom_sig_Banner %>% select(Uniq_ID)

## Merge nom sig with GW sig

Merged_1 <- merge(GW_sig_ROSMAP, Nom_sig_Banner, by="Uniq_ID")

## Test which SNPs are associated with three or fewer genes

SNPs_to_test <- Merged_1 %>% select(SNP)

SNPs_assoc_with_3_or_less <- SNPs_to_test %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT)) %>% filter(COUNT < 4)

## Final merge with ROSMAP and export

Final_merge <- merge(SNPs_assoc_with_3_or_less, Merged_1, by="SNP")

write.table(Final_merge, "~/Desktop/SZ_PES_mtCOJO_norm/MR_input/Protein/Brain/Processed_raw_input_ROSMAP_consistent_with_Banner_3_or_less_proteins.txt",
            sep = "\t", row.names = F, quote = F)


## Brain eQTLs - MetaBrain genome-wide significant

MetaBrain_raw <- fread("~/Desktop/SZ_PES_mtCOJO_norm/MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/All_chr_MetaBrain_GW_sig_eQTLs_Cortex.txt", header = T)

## Test which SNPs are associated with three or fewer genes

SNPs_to_test_eQTL <- MetaBrain_raw %>% select(SNP)

SNPs_assoc_with_3_or_less <- SNPs_to_test_eQTL %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT)) %>% filter(COUNT < 4)

## Merge and output

MetaBrain_output <- merge(MetaBrain_raw, SNPs_assoc_with_3_or_less, by="SNP")

write.table(MetaBrain_output, file="~/Desktop/SZ_PES_mtCOJO_norm/MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Less_than_4_genes_per_QTL_MetaBrain.txt", sep="\t", quote = F, row.names = F)

## Basal Ganglia

## Brain eQTLs - MetaBrain genome-wide significant

MetaBrain_raw_BG <- fread("~/Desktop/SZ_PES_mtCOJO_norm/MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/dbSNP_raw_GW_sig_basal_ganglia.txt", header = T)

## Test which SNPs are associated with three or fewer genes

SNPs_to_test_eQTL <- MetaBrain_raw_BG %>% select(rsid)

SNPs_assoc_with_3_or_less <- SNPs_to_test_eQTL %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT)) %>% filter(COUNT < 4)

MetaBrain_output_BG <- merge(MetaBrain_raw_BG, SNPs_assoc_with_3_or_less, by="rsid")

write.table(MetaBrain_output_BG, file="~/Desktop/SZ_PES_mtCOJO_norm/MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Basal_ganglia_Less_than_4_genes_per_QTL_MetaBrain.txt", sep="\t", quote = F, row.names = F)

## Cerebellum

## Brain eQTLs - MetaBrain genome-wide significant

MetaBrain_raw_CB <- fread("~/Desktop/SZ_PES_mtCOJO_norm/MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/dbSNP_raw_cerebellum_GW_sig.txt", header = T)

## Test which SNPs are associated with three or fewer genes

SNPs_to_test_eQTL <- MetaBrain_raw_CB %>% select(SNP)

SNPs_assoc_with_3_or_less <- SNPs_to_test_eQTL %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT)) %>% filter(COUNT < 4)

MetaBrain_output_CB <- merge(MetaBrain_raw_CB, SNPs_assoc_with_3_or_less, by="SNP")

write.table(MetaBrain_output_CB, file="~/Desktop/SZ_PES_mtCOJO_norm/MR_input/mRNA/Brain/MetaBrain_eQTL_SMR/Cerebellum_Less_than_4_genes_per_QTL_MetaBrain.txt", sep="\t", quote = F, row.names = F)

## Blood eQTLgen

Blood_raw <- fread("~/Desktop/SZ_PES_mtCOJO_norm/MR_input/mRNA/Blood/Cis_bloodeQTLgen_gw_sig_RAW.txt", header = T)

## Test which SNPs are associated with three or fewer genes

SNPs_to_test_eQTL <- Blood_raw %>% select(SNP)

SNPs_assoc_with_3_or_less <- SNPs_to_test_eQTL %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT)) %>% filter(COUNT < 4)

Blood_output <- merge(Blood_raw, SNPs_assoc_with_3_or_less, by="SNP")

write.table(Blood_output, file="~/Desktop/SZ_PES_mtCOJO_norm/MR_input/mRNA/Blood/RAW_no_LD_pruning_Blood_IVs_three_or_less_genes.txt", sep="\t", quote = F, row.names = F)
