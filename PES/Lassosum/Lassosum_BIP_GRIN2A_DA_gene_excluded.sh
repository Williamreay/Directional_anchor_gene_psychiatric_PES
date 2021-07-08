#!/bin/bash

#SBATCH --cpus-per-task=6
#SBATCH --time=130:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

while read GWAS_FILE SCORE_NAME;
  do
    lassosum --data $GWAS_FILE \
    --chr CHR --pos BP --A1 A1 --A2 A2 --pval P --n 51890 --OR OR \
    --ref ~/ukbb2/SZ_and_BIP_UKBB_subsets/SZ_cohort/Merged_SZ_UKBB_cohort --keep.ref ../PES/TRAINING_BIP_phenotypes/Kee_LDref_BIP.txt \
    --LDblocks EUR.hg19 --pheno ../PES/TRAINING_BIP_phenotypes/BIP_TRAINING_UKBB_pheno.txt \
    --covar ../PES/TRAINING_BIP_phenotypes/Covariates_TRAIN_BIP_UKBB.txt --test.bfile ~/ukbb2/SZ_and_BIP_UKBB_subsets/BIP_TRAIN_cohort/Merged_BIP_TRAIN_UKBB_cohort \
    --nthreads 6 --out BIP_TRAIN_lassosum_liberal/BIP_TRAIN_UKBB_LASSO_${SCORE_NAME};
done < DA_gene_removed/GRIN2A_list_BIP.txt

