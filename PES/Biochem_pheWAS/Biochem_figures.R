################################

## Effect of PES/PRS on blood/urine biochemical traits amongst UKBB participants who did not self-report mental illness in the MHQ

## Figures

## William Reay (2021)

#################################

library(ggplot2)
library(readxl)
library(dplyr)
library(data.table)
library(ggpubr)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem/")

## Volcano plot for FADS1 BIP and SZ

## SZ FADS1

SZ_FADS1 <- fread("Combined_results/SZ_SZ_FADS1_PES.txt", header = F)

SZ_FADS1 <- rename(SZ_FADS1, "field_ID"="V1", "Tstat"="V4", "P"="V5")

SZ_FADS1$threshold <- as.factor(ifelse(SZ_FADS1$P < 1.05e-04, 1, 0))

Biochem <- read_excel("Biochemical_traits_UKBB.xlsx")

SZ_FADS1 <- as.data.frame(merge(SZ_FADS1, Biochem, by = "field_ID"))

SZ_FADS1_plot <- ggplot(data = SZ_FADS1,
      aes(x=Tstat, y = -log10(P),
          colour = threshold, label = Biochem_trait)) +
  scale_color_manual(name="threshold", values = c("black","#3333FF")) +
  geom_point(alpha = 0.8, size = 1.75) +
  labs(x = expression("t value"), y = expression(paste("-log"[10], "P-value"))) +
  theme_bw() +
  geom_text(aes(label=ifelse(P < 1.05e-04 & Tstat < 0, as.character(Biochem_trait), '')), hjust=-0.1, vjust=-0.5) +
  geom_text(aes(label=ifelse(P < 1.05e-04 & Tstat > 0, as.character(Biochem_trait), '')), hjust=0.5, vjust=-0.5) +
  theme(legend.position ="none") +
  theme(legend.position ="none") +
  theme(axis.title = element_text(face="bold", size=12)) +
  ggtitle("SZ FADS1 network PES") +
  geom_hline(yintercept = 3.978811, linetype ="longdash") +
  xlim(c(-7,7)) +
  ylim(c(0,12))

## SZ PRS

SZ_PRS <- fread("Combined_results/SZ_SZ_PRS.txt", header = F)

SZ_PRS <- rename(SZ_PRS, "field_ID"="V1", "Tstat"="V4", "P"="V5")

SZ_PRS$threshold <- as.factor(ifelse(SZ_PRS$P < 1.05e-04, 1, 0))

SZ_PRS <- as.data.frame(merge(SZ_PRS, Biochem, by = "field_ID"))

PRS_plot <- ggplot(data = SZ_PRS,
                        aes(x=Tstat, y = -log10(P),
                            colour = threshold, label = Biochem_trait)) +
  scale_color_manual(name="threshold", values = c("black","orangered")) +
  geom_point(alpha = 0.8, size = 1.75) +
  labs(x = expression("t value"), y = expression(paste("-log"[10], " P-value"))) +
  theme_bw() +
  geom_text(aes(label=ifelse(P < 1.05e-04 & Tstat < 0, as.character(Biochem_trait), '')), hjust=-0.1, vjust=-0.5) +
  geom_text(aes(label=ifelse(P < 1.05e-04 & Tstat > 0, as.character(Biochem_trait), '')), hjust=0.5, vjust=-0.5) +
  theme(legend.position ="none") +
  theme(legend.position ="none") +
  theme(axis.title = element_text(face="bold", size=12)) +
  ggtitle("SZ PRS") +
  geom_hline(yintercept = 3.978811, linetype ="longdash") +
  xlim(c(-7,7)) +
  ylim(c(0,12))



