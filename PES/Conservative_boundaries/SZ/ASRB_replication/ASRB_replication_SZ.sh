#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=130:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

while read SCORE_FILE SCORE;
  do
    plink2 --bfile /home/control/data/genomic_data/arrays/610k/clean_rm_mhc/610K_MHC_removed \
    --score $SCORE_FILE --out ASRB_replication_full_cohort_ex_ASRB/ASRB_replication_SZ_${SCORE};
  done < ASRB_replication_input.txt
