#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=13:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## BIP - brain

for chr in $(seq 1 22);
do
Rscript ~/data/users/william/fusion_twas-master/FUSION.assoc_test.R \
--sumstats ~/data/users/william/Trans_diagnostic_repurposing_2020/Munged_dbSNP_sumstats/Munged_BIP_2019.sumstats.gz \
--weights  /home/control/data/users/william/PWAS_weights/Brain/ROSMAP_and_Banner_PWAS_weights/ROSMAP.n376.fusion.WEIGHTS/train_weights.pos \
--weights_dir /home/control/data/users/william/PWAS_weights/Brain/ROSMAP_and_Banner_PWAS_weights/ROSMAP.n376.fusion.WEIGHTS/ \
--ref_ld_chr ~/data/users/william/fusion_twas-master/LDREF/1000G.EUR. \
--chr ${chr} \
--out ROSMAP_brain_BIP/BIP_ROSMAP_PWAS_chr_${chr}
done
