#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=130:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## SZ UKBB - PES - traditional PES

while read GWAS_FILE SCORE_NAME THRESHOLD;
do
  python ../../PES/PES_generation_no_stats.py \
  --pathway_gwasfile $GWAS_FILE \
  --bfile ~/ukbb2/SZ_and_BIP_UKBB_subsets/SZ_cohort/Merged_SZ_UKBB_cohort \
  --PRsice2_binary PRSice_linux \
  --phenotype SZ_UKBB \
  --pathway $SCORE_NAME \
  --A1 A1 --A2 A2 --CHR CHR --BP BP --SNP SNP --P P \
  --p_threshold $THRESHOLD \
  --statistic OR \
  --target_type T;
  done < SZ_UKBB_input.txt
