########################################

## eQTL and pQTL MR

## William Reay (2021)

#######################################

library(TwoSampleMR)
library(dplyr)
library(data.table)


Harmonise_and_MR <- function(IVs, Outcome, Phenotype) {
  Extract <- format_data(dat = Outcome, type = "outcome",
                         snps = IVs$SNP, snp_col = "SNP", beta_col = "Beta",
                         se_col = "SE", effect_allele_col = "A1",
                         other_allele_col = "A2", pval_col = "P")
  Harmonised_phenotype <- harmonise_data(IVs, Extract, action = 3)
  MR <- mr(Harmonised_phenotype, method_list = c("mr_ivw_fe", "mr_wald_ratio"))
  return(MR)
}

save(Harmonise_and_MR, file="~/Desktop/SZ_PES_mtCOJO_norm/MR_input/MR_func.Rda")
