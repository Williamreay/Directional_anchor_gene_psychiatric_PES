#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=12:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

#Schizophrenia BIP GTExv7 brain

while read PHENOTYPE TISSUE CHROMOSOME DIR;
do
Rscript ~/data/users/william/fusion_twas-master/FUSION.assoc_test.R \
--sumstats  ~/data/users/william/Trans_diagnostic_repurposing_2020/Munged_dbSNP_sumstats/Munged_BIP_2019.sumstats.gz \
--weights /home/control/data/users/william/GTEx_TWAS_weights/${DIR}/${TISSUE}.P01.pos \
--weights_dir /home/control/data/users/william/GTEx_TWAS_weights/${DIR} \
--ref_ld_chr ~/data/users/william/fusion_twas-master/LDREF/1000G.EUR. \
--chr $CHROMOSOME --out ~/data/users/william/Trans_diagnostic_repurposing_2020/TWAS/BIP_2019/${PHENOTYPE}_${TISSUE}_${CHROMOSOME};
done < BIP_input.txt
