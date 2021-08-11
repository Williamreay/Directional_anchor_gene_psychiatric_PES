#############################

## Correlation of PES/PRS with self-reported MHQ phenotypes

## SZ and BIP excluded + controls in the discovery analyses

## Sensitivity analyses

## William Reay (2021)

############################

library(dplyr)
library(data.table)
library(easyGgplot2)
library(corrplot)
library(reshape2)
library(ComplexHeatmap)
library(circlize)


setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_self_reported_mental_illness/")

## Read in PES/PRS

PES_PRS <- fread("../UKBB_MHQ_no_self_reported_biochem/Entire_EUR_genotyped_UKBB_scores/FINAL_SCORES_PES_PRS/All_SZ_BIP_PES_PRS_all_UKBB.txt", header = T)

PES_PRS <- PES_PRS[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28)]

## Load in MHQ self-report data

MHQ_self_report <- fread("UKBB_self_report_mental_illness_df_SZ_BIP_excl.txt", header = T)

MHQ_self_report <- rename(MHQ_self_report, "IID"="f.eid")

## Merge

Merged_MHQ_scoring <- merge(PES_PRS, MHQ_self_report, by = "IID")

Merged_MHQ_scoring$Batch <- as.factor(Merged_MHQ_scoring$Batch)

## Define list of scores

Scores_to_test <- as.list(colnames(Merged_MHQ_scoring[,2:15]))

## SZ CACNA1C PES and depression

Dep <- Merged_MHQ_scoring %>% filter(No_mental_health == 0 | Depression == 1)

Dep$SZ_CACNA1C_PES <- as.numeric(scale(Dep$SZ_CACNA1C_PES))
Dep$SZ_PRS <- as.numeric(scale(Dep$SZ_PRS))

Dep_test <- glm(Depression ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
  PC8 + PC9 + PC10 + Batch + SZ_PRS + SZ_CACNA1C_PES, family = "binomial",
  data = Dep)

anova(Dep_test, test="Chisq")

## SZ RPS17 PES and OCD

OCD <- Merged_MHQ_scoring %>% filter(No_mental_health == 0 | OCD == 1)

OCD$SZ_RPS17_PES <- as.numeric(scale(OCD$SZ_RPS17_PES))
OCD$SZ_PRS <- as.numeric(scale(OCD$SZ_PRS))

OCD_test <- glm(OCD ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                  PC8 + PC9 + PC10 + Batch + SZ_PRS + SZ_RPS17_PES, family = "binomial",
                data = OCD)

anova(OCD_test, test="Chisq")

## Make heatmap of self-report Z scores

Hm_input <- fread("FDR_FWER_MHQ_phenotypes.csv", header = T)

Hm_input <- Hm_input %>% select(z, Phenotype, Score)

Input <- dcast(Hm_input, Phenotype ~ Score, value.var = "z")

Input_mat <- as.matrix(Input)

rownames(Input_mat) <- c("ADHD", "Agoraphobia", "Anorexia nervosa",
                         "Anxiety/GAD", "ASD", "Bulimia nervosa", 
                         "Depression", "OCD", "Other phobia", "Other psychosis",
                         "Over eating/binge eating", "Paniac attack", "Personality disorder",
                         "Social anxiety/Social phobia")

colnames(Input_mat) <- c("BIP CACNA1C PES", "BIP FADS1 PES",
                         "BIP FES PES", "BIP GRIN2A PES", "BIP PCCB PES",
                         "BIP PRS", "BIP RPS17 PES", "SZ CACNA1C PES",
                         "SZ FADS1 PES", "SZ FES PES", "SZ GRIN2A PES",
                         "SZ PCCB PES", "SZ PRS", "SZ RPS17 PES")
                        

Input_mat <- Input_mat[, -1]

Input_mat <- `dimnames<-`(`dim<-`(as.numeric(Input_mat), dim(Input_mat)), dimnames(Input_mat))

Heatmap(Input_mat, name = "Z", rect_gp = gpar(col = "white", lwd = 2), 
        show_column_dend = TRUE, show_row_dend = TRUE, clustering_distance_rows = "pearson",
        row_dend_width = unit(1, "cm"),
        row_names_gp = grid::gpar(fontsize = 8.5), column_names_gp = grid::gpar(fontsize = 7.5),
        column_names_rot = 50, column_title_gp = grid::gpar(fontsize = 11.5, fontface="bold"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(Input_mat[i, j] > 3.657) {
            grid.text("***", x, y)
          }
          if(Input_mat[i, j] > 2.77 && Input_mat[i, j] < 3.657) {
            grid.text("**", x, y)
          }
          if(Input_mat[i, j] > 1.96 && Input_mat[i, j] < 2.77) {
            grid.text("*", x, y)
          }
          if(Input_mat[i, j] < -1.96) {
            grid.text("*", x, y)
          }
        }
        )









