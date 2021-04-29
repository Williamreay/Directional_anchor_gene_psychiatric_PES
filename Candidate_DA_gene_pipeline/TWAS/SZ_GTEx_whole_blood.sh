#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=12:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

#Schizophrenia  GTEx v7 whole blood

while read PHENOTYPE TISSUE CHROMOSOME DIR;
do
Rscript ~/data/users/william/fusion_twas-master/FUSION.assoc_test.R \
--sumstats  ~/data/users/william/Trans_diagnostic_repurposing_2020/Munged_dbSNP_sumstats/Munged_SZ_PGC3.sumstats.gz \
--weights /home/control/data/users/william/GTEx_TWAS_weights/Whole_Blood/Whole_Blood.P01.pos \
--weights_dir /home/control/data/users/william/GTEx_TWAS_weights/ \
--ref_ld_chr ~/data/users/william/fusion_twas-master/LDREF/1000G.EUR. \
--chr $CHROMOSOME --out ~/data/users/william/Trans_diagnostic_repurposing_2020/TWAS/SZ_PGC3_2020/${PHENOTYPE}_${TISSUE}_${CHROMOSOME};
done < SZ_input_2.txt
