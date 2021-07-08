##########################################

## STRING networks for candidate DA-genes - investigation of random walk algorithim

## William Reay (2021)

##########################################

set.seed(676742)

library(dplyr)
library(data.table)


setwd("~/Desktop/SZ_PES_mtCOJO_norm/Candidate_DAG/")

## Read in list of candidate DA-genes

Candidates <- fread("Final_list_of_candidate_DAG_SZ_and_BIP.txt", header = F)

## Load STRING database - version 11

STRINGdb_v11 <- fread("9606.protein.links.detailed.v11.0.txt.gz", header = T)

## Filter high confidence edges (experimental or database > 0.7) - based on experimental and database evidence

High_conf_STRINGdb_v11_exp_db <- STRINGdb_v11 %>% select(protein1, protein2, experimental, database) %>% filter(experimental > 700 | database > 700) 

## Read in list of database identifiers for STRING

db_identifiers <- fread("9606.protein.info.v11.0.txt.gz", header = T)

db_identifiers <- db_identifiers %>% select(protein_external_id, preferred_name)

## Map the input DA-genes to their STRING identifier

Candidates <- rename(Candidates, "preferred_name"="V1")

Merged_Candidates <- merge(Candidates, db_identifiers, by = "preferred_name")

## Define function to extract candidates

STRING_extract <- function(db, gene) {
  db_extract <- db %>% filter(protein1 == gene)
  return(db_extract)
  
}

FES <- STRING_extract(High_conf_STRINGdb_v11_exp_db, "9606.ENSP00000331504")

## Define list of seed genes and extract from high confidence network

List_seed_genes <- as.list(Merged_Candidates$protein_external_id)

Output <- lapply(List_seed_genes, STRING_extract, db = High_conf_STRINGdb_v11_exp_db)

## Rename the STRING IDs protein IDs

CACNA1C <- Output[[1]]

FADS1 <- Output[[2]]

FES <- Output[[3]]

GRIN2A <- Output[[4]]

PCCB <- Output[[5]]

RPS17 <- Output[[6]]

## Combine and rename

Merge_and_rename <- function(gene, df) {
  
  df_rename <- rbind(df[1,1], df$protein2, use.names = F)
  
  df_rename <- rename(df_rename, "protein_external_id"="protein1")
  
  gene_merged <- merge(df_rename, db_identifiers, by="protein_external_id")
}

FINAL_CACNA1C <- Merge_and_rename(CACNA1C, CACNA1C)

write.table(FINAL_CACNA1C, file="STRING_networks/CACNA1C_network.txt",
            sep = "\t", row.names = F, quote = F)

FINAL_FADS1 <- Merge_and_rename(FADS1, FADS1)

write.table(FINAL_FADS1, file="STRING_networks/FADS1_network.txt",
            sep = "\t", row.names = F, quote = F)

FINAL_GRIN2A <- Merge_and_rename(GRIN2A, GRIN2A)

write.table(FINAL_GRIN2A, file="STRING_networks/GRIN2A_network.txt",
            sep = "\t", row.names = F, quote = F)

FINAL_FES <- Merge_and_rename(FES, FES)

write.table(FINAL_FES, file="STRING_networks/FES_network.txt",
            sep = "\t", row.names = F, quote = F)

FINAL_PCCB <- Merge_and_rename(PCCB, PCCB)

write.table(FINAL_PCCB, file="STRING_networks/PCCB_network.txt",
            sep = "\t", row.names = F, quote = F)

FINAL_RPS17 <- Merge_and_rename(RPS17, RPS17)

write.table(FINAL_RPS17, file="STRING_networks/RPS17_network.txt",
            sep = "\t", row.names = F, quote = F)

