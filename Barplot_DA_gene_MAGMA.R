########################

## Bar plots of DA-gene network associations

## William Reay (2021)

########################

library(readxl)
library(ggplot2)
library(ggpubr)

## Conservative genic boundaires

Cons <- read_excel("~/Desktop/SZ_PES_mtCOJO_norm/Updated_SZ_BIP_PES_pathways/DA_gene_networks_GSA/DA_gene_network_assoc_FP_CONSERVATIVE.xlsx")

C <- ggplot(data = Cons, aes(x=Set, y=-log10(P), fill = as.factor(Gene_included))) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  facet_wrap(~Phenotype,strip.position="top",nrow=3,scales = "free_y") +
  coord_flip() + theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(size=10)) +
  geom_hline(yintercept=1.30103, lty=2) +
  ylab("-log10 MAGMA P-value") + xlab(" ") +
  ggtitle("Conservative genic boundaries") +
  ylim(0,9)

  
Liberal <- read_excel("~/Desktop/SZ_PES_mtCOJO_norm/Updated_SZ_BIP_PES_pathways/DA_gene_networks_GSA/DA_gene_networks_LIBERAL.xlsx")

L <- ggplot(data = Liberal, aes(x=Set, y=-log10(P), fill = Gene_included)) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  facet_wrap(~Phenotype,strip.position="top",nrow=3,scales = "free_y") +
  coord_flip() + theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(size=10)) +
  geom_hline(yintercept=1.30103, lty=2) +
  ylab("-log10 MAGMA P-value") + xlab(" ") +
  ggtitle("Liberal genic boundaries") +
  ylim(0,9)

All <- ggarrange(C, L, nrow = 1, ncol=2, common.legend = T, legend = "bottom")
  
 