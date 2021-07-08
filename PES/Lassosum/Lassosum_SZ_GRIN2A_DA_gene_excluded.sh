#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=130:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

while read GWAS_FILE SCORE_NAME;
  do
    lassosum --data $GWAS_FILE \
    --chr CHR --pos BP --A1 A1 --A2 A2 --pval P --n 161405 --OR OR \
    --ref ~/ukbb2/SZ_and_BIP_UKBB_subsets/BIP_TRAIN_cohort/Merged_BIP_TRAIN_UKBB_cohort \
    --LDblocks EUR.hg19 --pheno ../SZ_phenotypes/SZ_UKBB_pheno.txt --keep.ref ../SZ_phenotypes/SZ_lassosum_ref_keep.txt \
    --covar ../SZ_phenotypes/Covariates_SZ_UKBB.txt --test.bfile ~/ukbb2/SZ_and_BIP_UKBB_subsets/SZ_cohort/Merged_SZ_UKBB_cohort \
    --nthreads 6 --out SZ_UKBB_lassosum_liberal/SZ_UKBB_LASSO_${SCORE_NAME};
done < DA_gene_removed/GRIN2A_list.txt
