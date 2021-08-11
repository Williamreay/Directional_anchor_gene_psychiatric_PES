######################################

## Z score test of sexual dimorphism between PES effects on the biochemical traits

## William Reay (2021)

######################################

library(data.table)
library(dplyr)
library(ggplot2)
library(readxl)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem/Sex_stratified/")

## Read in male results

Male <- read_excel("All_male_combined.xlsx")

## Read in female results

Female <- read_excel("All_female_combined.xlsx")

## Combine df

Merged_both_sexes <- merge(Male, Female, by=c("field_ID", "Score"))

## Perform Z test - c.f. https://doi.org/10.1016/j.biopsych.2020.12.024

Merged_both_sexes$Z_dimorphism <- ((Merged_both_sexes$Female_beta - Merged_both_sexes$Male_beta)/sqrt(Merged_both_sexes$Female_se^2+Merged_both_sexes$Male_se^2))

Merged_both_sexes$P_dimorphism <- pnorm(abs(Merged_both_sexes$Z_dimorphism), lower.tail=FALSE) * 2

## Identify P values < 0.05

Significant_dimorphism <- Merged_both_sexes %>% filter(P_dimorphism < 0.05)

write.csv(Significant_dimorphism, file="../Biochem_sensitivity_analyses/Nominal_sexual_dimorphism_evidence.csv",
          row.names = F, quote = F)

Biochem <- read_excel("../Biochemical_traits_UKBB.xlsx")

Significant_dimorphism <- rename(Significant_dimorphism, "field_ID"="Biochem")

Significant_dimorphism <- merge(Significant_dimorphism, Biochem, by="field_ID")

Significant_dimorphism$Abs_Z <- abs(Significant_dimorphism$Z_dimorphism)

## Order and retain largest value

Significant_dimorphism <- as.data.frame(Significant_dimorphism %>% group_by(field_ID) %>% top_n(1, Abs_Z))

## Plot Z

ggplot(data = Significant_dimorphism, aes(x=field_ID, y=Z_dimorphism, fill = field_ID)) +
  geom_bar(stat="identity", color="black") +
  coord_flip() +
  theme_bw() +
  ylab("Z score of sex difference") +
  xlab(" ") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype ="longdash")

