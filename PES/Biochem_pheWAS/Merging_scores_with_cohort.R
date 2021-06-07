################################

## Merging scores with UKBB participants who did not self-report mental illness in the MHQ

## SZ and BIP cohorts excluded

## William Reay (2021)

#################################

library(dplyr)
library(data.table)

## Load in merged scores

All_scores <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem/Entire_EUR_genotyped_UKBB_scores/FINAL_SCORES_PES_PRS/All_SZ_BIP_PES_PRS_all_UKBB.txt",
                    sep ="\t")

Colnames_all_scores <- make.unique(names(All_scores))

colnames(All_scores) <- Colnames_all_scores

## Load in MHQ participants with no self-report not in SZ or BIP

MHQ_IDs <- fread("~/Desktop/SZ_PES_mtCOJO_norm/UKBB_phenotypes/UKBB_biochemical_pheWAS_IDs/Measured_MHQ_controls_biochem.txt")

MHQ_IDs <- rename(MHQ_IDs, "IID"="f.eid")

## Merge

Final_merged <- merge(MHQ_IDs, All_scores, by ="IID")

write.table(Final_merged,
            file ="~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem/Merged_all_scores_all_biochem.txt",
            sep = "\t", row.names = F, quote = F)

