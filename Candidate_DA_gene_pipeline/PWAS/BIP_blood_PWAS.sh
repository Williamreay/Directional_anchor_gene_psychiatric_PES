#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=18:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## Bipolar blood PWAS 

for chr in $(seq 1 22);
do
Rscript Blood/scripts/PWAS.assoc_test.R \
--sumstats ~/data/users/william/Trans_diagnostic_repurposing_2020/Munged_dbSNP_sumstats/Munged_BIP_2019.sumstats.gz \
--weights Blood/PWAS_EA/Plasma_Protein_EA_hg19.pos \
--weights_dir Blood/PWAS_EA/Plasma_Protein_weights_EA_enet \
--force_model enet \
--ref_ld_chr Blood/LDref/EUR/chr \
--chr $chr \
--out Blood_BIP_results/BIP_PWAS_enet_chr_chr${chr};
done
