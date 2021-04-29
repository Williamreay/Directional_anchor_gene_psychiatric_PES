#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=19:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## BIP 2019 FUSION TWAS

## PsychENCODE TWAS weights

## William Reay

for chr in $(seq 1 22);
do
Rscript /home/control/data/users/william/fusion_twas-master/FUSION.assoc_test.R \
--sumstats ~/data/users/william/Trans_diagnostic_repurposing_2020/Munged_CHR_BP_sumstats/Munged_CHR_BP_BIP_2019.sumstats.gz \
--weights /home/control/data/users/william/PEC_TWAS/PEC_WEIGHTS.pos \
--weights_dir /home/control/data/users/william/PEC_TWAS/ \
--ref_ld_chr /home/control/data/users/william/CHR_BP_g1000_LDREF/All_common_SNPs_1000G.EUR. \
--chr $chr --out /home/control/data/users/william/Trans_diagnostic_repurposing_2020/TWAS/BIP_2019/BIP_2019_chr${chr}_PEC.txt;
done
