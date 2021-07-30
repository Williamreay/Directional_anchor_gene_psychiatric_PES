#########################################

## Concatenating chromosome-wise PES/PRS in the entire EUR UKBB

## William Reay (2021)

#########################################

library(data.table)
library(dplyr)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem/Entire_EUR_genotyped_UKBB_scores/")

## Define concatenation function

PES_concat <- function(df, name){
  ## Format df with unique colnames
  Colnames_df <- make.unique(names(df))
  colnames(df) <- Colnames_df
  ## Select columns for PES or PRS
  Score_df <- df %>% select(IID, contains("AVG"))
  ## Sum chromosome-wise score
  Score_df[[name]] <- rowSums(Score_df[,-1])
  Score_df <- Score_df %>% select(IID, name)
  return(Score_df)
}

#### SZ ####

## SZ PRS

SZ_PRS_raw <- fread("UKBB_SZ_PRS", header = T, sep="\t")

SZ_PRS_concat <- PES_concat(SZ_PRS_raw, "SZ_PRS")

write.table(SZ_PRS_concat, file="FINAL_SCORES_PES_PRS/SZ_PRS_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## CACNA1C SZ PES

CACNA1C_PES_raw <- fread("UKBB_SZ_CACNA1C", header = T, sep="\t")

SZ_CACNA1C_concat <- PES_concat(CACNA1C_PES_raw, "SZ_CACNA1C_PES")

write.table(SZ_CACNA1C_concat, file="FINAL_SCORES_PES_PRS/SZ_CACNA1C_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## GRIN2A SZ PES

GRIN2A_PES_raw <- fread("UKBB_SZ_GRIN2A", header = T, sep="\t")

SZ_GRIN2A_concat <- PES_concat(GRIN2A_PES_raw, "SZ_GRIN2A_PES")

write.table(SZ_GRIN2A_concat, file="FINAL_SCORES_PES_PRS/SZ_GRIN2A_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## FES SZ PES

FES_PES_raw <- fread("UKBB_SZ_FES", header = T, sep="\t")

SZ_FES_concat <- PES_concat(FES_PES_raw, "SZ_FES_PES")

write.table(SZ_FES_concat, file="FINAL_SCORES_PES_PRS/SZ_FES_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## FADS1 SZ PES

FADS1_PES_raw <- fread("UKBB_SZ_FADS1", header = T, sep="\t")

SZ_FADS1_concat <- PES_concat(FADS1_PES_raw, "SZ_FADS1_PES")

write.table(SZ_FADS1_concat, file="FINAL_SCORES_PES_PRS/SZ_FADS1_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## PCCB SZ PES

PCCB_PES_raw <- fread("UKBB_SZ_PCCB", header = T, sep="\t")

SZ_PCCB_concat <- PES_concat(PCCB_PES_raw, "SZ_PCCB_PES")

write.table(SZ_PCCB_concat, file="FINAL_SCORES_PES_PRS/SZ_PCCB_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## RPS17 SZ PES

RPS17_PES_raw <- fread("UKBB_SZ_RPS17", header = T, sep="\t")

SZ_RPS17_concat <- PES_concat(RPS17_PES_raw, "SZ_RPS17_PES")

write.table(SZ_RPS17_concat, file="FINAL_SCORES_PES_PRS/SZ_RPS17_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

##### BIP #####

BIP_PRS_raw <- fread("UKBB_BIP_PRS", header = T, sep="\t")

BIP_PRS_concat <- PES_concat(BIP_PRS_raw, "BIP_PRS")

write.table(BIP_PRS_concat, file="FINAL_SCORES_PES_PRS/BIP_PRS_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## CACNA1C BIP PES

CACNA1C_PES_raw <- fread("UKBB_BIP_CACNA1C", header = T, sep="\t")

BIP_CACNA1C_concat <- PES_concat(CACNA1C_PES_raw, "BIP_CACNA1C_PES")

write.table(BIP_CACNA1C_concat, file="FINAL_SCORES_PES_PRS/BIP_CACNA1C_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## GRIN2A BIP PES

GRIN2A_PES_raw <- fread("UKBB_BIP_GRIN2A", header = T, sep="\t")

BIP_GRIN2A_concat <- PES_concat(GRIN2A_PES_raw, "BIP_GRIN2A_PES")

write.table(BIP_GRIN2A_concat, file="FINAL_SCORES_PES_PRS/BIP_GRIN2A_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## FES BIP PES

FES_PES_raw <- fread("UKBB_BIP_FES", header = T, sep="\t")

BIP_FES_concat <- PES_concat(FES_PES_raw, "BIP_FES_PES")

write.table(BIP_FES_concat, file="FINAL_SCORES_PES_PRS/BIP_FES_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## FADS1 BIP PES

FADS1_PES_raw <- fread("UKBB_BIP_FADS1", header = T, sep="\t")

BIP_FADS1_concat <- PES_concat(FADS1_PES_raw, "BIP_FADS1_PES")

write.table(BIP_FADS1_concat, file="FINAL_SCORES_PES_PRS/BIP_FADS1_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## PCCB BIP PES

PCCB_PES_raw <- fread("UKBB_BIP_PCCB", header = T, sep="\t")

BIP_PCCB_concat <- PES_concat(PCCB_PES_raw, "BIP_PCCB_PES")

write.table(BIP_PCCB_concat, file="FINAL_SCORES_PES_PRS/BIP_PCCB_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)

## RPS17 BIP PES

RPS17_PES_raw <- fread("UKBB_BIP_RPS17", header = T, sep="\t")

BIP_RPS17_concat <- PES_concat(RPS17_PES_raw, "BIP_RPS17_PES")

write.table(BIP_RPS17_concat, file="FINAL_SCORES_PES_PRS/BIP_RPS17_PES_UKBB.txt",
            sep = "\t", row.names = F, quote = F)
