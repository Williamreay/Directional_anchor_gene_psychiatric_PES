
library(dplyr)
library(data.table)
library(readxl)

setwd("/Users/williamreay/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem")

All <- read_excel("Combined_results/All_biochem_combined_results.xlsx")

Biochem_names <- fread("Biochem_names.txt",header = T)

Merged <- merge(All, Biochem_names, by="field_ID")

Merged$FDR <- p.adjust(Merged$P, method="fdr")
Merged$FWER <- p.adjust(Merged$P, method="bonferroni")

write.table(Merged, file="Final_FWER_FDR_biochem_pheWAS_both_sexes.txt",
            sep = "\t",row.names = F, quote = F)

## Female stratified

Female <- read_excel("Sex_stratified/All_female_combined.xlsx")

Merged_fe <- merge(Female, Biochem_names, by="field_ID")

Merged_fe$FDR <- p.adjust(Merged_fe$P, method="fdr")
Merged_fe$FWER <- p.adjust(Merged_fe$P, method="bonferroni")

write.table(Merged_fe, file="Sex_stratified/Final_FWER_FDR_biochem_pheWAS_FEMALES.txt",
            sep = "\t",row.names = F, quote = F)

## Male stratified

Male <- read_excel("Sex_stratified/All_male_combined.xlsx")

Merged_male <- merge(Male, Biochem_names, by="field_ID")

Merged_male$FDR <- p.adjust(Merged_male$P, method="fdr")
Merged_male$FWER <- p.adjust(Merged_male$P, method="bonferroni")

write.table(Merged_male, file="Sex_stratified/Final_FWER_FDR_biochem_pheWAS_MALES.txt",
            sep = "\t",row.names = F, quote = F)

