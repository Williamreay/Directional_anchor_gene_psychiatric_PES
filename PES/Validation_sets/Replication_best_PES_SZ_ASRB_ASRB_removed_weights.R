#############################################

## Association of PES and PRS with SZ

## Replication in the ASRB using best performing summary statistics

## Best performing PES

## William Reay (2021)

#############################################

setwd("~/Desktop/SZ_PES_mtCOJO_norm/")

library(data.table)
library(dplyr)
library(rcompanion)
library(corrplot)
library(ggplot2)
library(readxl)
library(viridis)

## Read in concatenated PES/PRS file

All_SZ_scores <- fread("Best_DA_gene_PES/SZ_best_con_or_lib/ASRB_SZ_VALIDATION/Combined_ASRB_replication.txt", header = T, sep = "\t")

Colnames_all_SZ_scores <- make.unique(names(All_SZ_scores))

colnames(All_SZ_scores) <- Colnames_all_SZ_scores

## Read in phenotype and covariate data

SZ_pheno <- read.csv("Best_DA_gene_PES/SZ_best_con_or_lib/ASRB_SZ_VALIDATION/Merged_ASRB_case_control.csv", header = T)

## Merge all three df

SZ_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                     list(All_SZ_scores, SZ_pheno))

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

SZ_merged[,c(6, 12, 18, 24, 30, 36, 42)] <- lapply(SZ_merged[,c(6, 12, 18, 24, 30, 36, 42)], function(x) c(scale(x)))

## Define function to test each score by itself and extract beta, SE, and P

## List of sets to test

Scores_to_test <- as.list(colnames(SZ_merged[,c(6, 12, 18, 24, 30, 36, 42)]))

## Function

PES_PRS_SZ_ASRB <- function(v){
  Score_model <- glm(glue::glue('PHENOTYPE ~ SEX + age + PC1 + PC2 + PC3 + {v}'), family = "binomial", data = SZ_merged)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

SZ_score <- sapply(Scores_to_test, PES_PRS_SZ_ASRB)

## Extract beta, se, z, and p value for each gene

SZ_extract <- apply(SZ_score, 2, function(x) return(as.data.frame(x$coefficients)[7,1:4]))

SZ_results <- data.frame()

for (i in 1:length(SZ_extract)) {
  SZ_results <- rbind(SZ_results, SZ_extract[[i]])
}
rownames(SZ_results) <- Scores_to_test

write.csv(SZ_results, file="Best_DA_gene_PES/SZ_best_con_or_lib/ASRB_SZ_VALIDATION/Results/VALIDATION_SZ_ASRB_best_PES_PRS_results.csv", quote = F, row.names = T)


## Make forest plot of SZ OR per SD in the training and validation set

SZ_results$Score <- Scores_to_test

SZ_results$Cohort <- "Validation"

## Read in training results for best PES

SZ_train <- read_excel("Best_DA_gene_PES/SZ_best_PES_TRAIN_summary.xlsx")

SZ_train <- SZ_train %>% select(Score, Beta, SE)

SZ_train$Cohort <- "Training"

SZ_results <- SZ_results %>% select(Score, Estimate, 'Std. Error', Cohort)

SZ_results <- rename(SZ_results, "Beta"="Estimate", "SE"="Std. Error")

Plot_df <- rbind(SZ_train, SZ_results)


Plot_df$Score_name <- c("PRS", "CACNA1C netowrk PES", "FADS1 network PES", 
                        "FES network PES", "GRIN2A network PES",
                        "PCCB network PES", "RPS17 network PES", "CACNA1C netowrk PES", "FADS1 network PES", 
                        "FES network PES", "GRIN2A network PES",
                        "PCCB network PES", "PRS", "RPS17 network PES")

## Plot results

Plot_df$OR <- exp(Plot_df$Beta)
Plot_df$UOR <- Plot_df$OR + 1.96*Plot_df$SE
Plot_df$LOR <- Plot_df$OR - 1.96*Plot_df$SE

Plot_df <- Plot_df %>% filter(Score_name !="PRS")

FP_SZ <- ggplot(data = Plot_df, aes(x=Score_name, y=OR, ymin=LOR, ymax=UOR, colour=Score_name)) +
  geom_pointrange() +
  geom_hline(yintercept=1, lty=2) +
  coord_flip() +
  ylab("Schizophrenia OR [95% CI] per SD increase in PES") +
  theme_bw() +
  theme(legend.position = "null", axis.title.y = element_blank()) +
  facet_wrap(~Cohort,strip.position="top",nrow=2,scales = "free_y")

FP_SZ + scale_color_brewer(palette = "Paired")

## Test association of PES/PRS with cognitive deficit relative to cognitively spared

CD_CS_df <- read.csv("Best_DA_gene_PES/SZ_best_con_or_lib/ASRB_SZ_VALIDATION/CD_vs_CS_GoM1.csv", header = T)
  
  
SZ_merged_GoM <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                      list(All_SZ_scores, CD_CS_df, SZ_pheno))

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

SZ_merged_GoM[,c(6, 12, 18, 24, 30, 36, 42)] <- lapply(SZ_merged_GoM[,c(6, 12, 18, 24, 30, 36, 42)], function(x) c(scale(x)))

## Function

PES_PRS_SZ_ASRB_GoM <- function(v){
  Score_model <- glm(glue::glue('GoM_coded ~ SEX.x + age + PC1.x + PC2.x + PC3.x + {v}'), family = "binomial", data = SZ_merged_GoM)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

SZ_score <- sapply(Scores_to_test, PES_PRS_SZ_ASRB_GoM)

## Extract beta, se, z, and p value for each gene

SZ_extract <- apply(SZ_score, 2, function(x) return(as.data.frame(x$coefficients)[7,1:4]))

SZ_results <- data.frame()

for (i in 1:length(SZ_extract)) {
  SZ_results <- rbind(SZ_results, SZ_extract[[i]])
}
rownames(SZ_results) <- Scores_to_test

write.csv(SZ_results, file="Best_DA_gene_PES/SZ_best_con_or_lib/ASRB_SZ_VALIDATION/Results/GoM_CD_CS_ASRB_best_PES_PRS_results.csv", quote = F, row.names = T)

## Test association with clozapine prescription - proxy for treatment resistance

Cloz_df <- fread("Best_DA_gene_PES/SZ_best_con_or_lib/ASRB_SZ_VALIDATION/Clozapine_ASRB_new_df.txt", header = T)

SZ_merged_Cloz <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                        list(All_SZ_scores, Cloz_df, SZ_pheno))

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

SZ_merged_Cloz[,c(6, 12, 18, 24, 30, 36, 42)] <- lapply(SZ_merged_Cloz[,c(6, 12, 18, 24, 30, 36, 42)], function(x) c(scale(x)))

## Function

PES_PRS_SZ_ASRB_Cloz <- function(v){
  Score_model <- glm(glue::glue('Clozapine ~ SEX + age + PC1 + PC2 + PC3 + {v}'), family = "binomial", data = SZ_merged_Cloz)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

SZ_score <- sapply(Scores_to_test, PES_PRS_SZ_ASRB_Cloz)

## Extract beta, se, z, and p value for each gene

SZ_extract <- apply(SZ_score, 2, function(x) return(as.data.frame(x$coefficients)[7,1:4]))

SZ_results <- data.frame()

for (i in 1:length(SZ_extract)) {
  SZ_results <- rbind(SZ_results, SZ_extract[[i]])
}
rownames(SZ_results) <- Scores_to_test

write.csv(SZ_results, file="Best_DA_gene_PES/SZ_best_con_or_lib/ASRB_SZ_VALIDATION/Results/Clozapine_ASRB_best_PES_PRS_results.csv", quote = F, row.names = T)

## Identify number of cases with top decile PES

Decile_df <- SZ_merged %>%
  mutate_at(vars(ends_with("AVG")), .funs = list(~ntile(.,10)))

## Assign individuals in the top PES decile

Top_decile_PES <- Decile_df %>%
  mutate_at(vars(ends_with("AVG")), .funs = list(~if_else(. == 10,1, 0)))

Top_decile_PES$PES_decile_sum <- rowSums(Top_decile_PES[, c(6, 12, 18, 24, 30, 42)])

Top_decile_PES$Binary_decile <- ifelse(Top_decile_PES$PES_decile_sum > 0, 1, 0)

table(Top_decile_PES$SZ, Top_decile_PES$PES_decile_sum)

Dec_tab <- table(Top_decile_PES$SZ, Top_decile_PES$Binary_decile)

Decile_SZ <- glm(PHENOTYPE ~  Binary_decile + SZ_PRS_AVG + SZ_PRS_AVG + PC1 + PC2 + PC3 + SEX + age,
                  family = "binomial", data = Top_decile_PES)

PRS_decile <- glm(Binary_decile ~ PRS,
                  family = "binomial", data = Top_decile_PES)

## Make mosaic plot of SZ cases and controls per decile grouping

mosaic( ~ SZ + Binary_decile, data = Top_decile_PES,
        highlighting = "SZ", highlighting_fill = c("grey", "plum"))  

Bottom_decile_PRS <- Decile_df %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))

Bottom_decile_PRS_only <- Bottom_decile_PRS %>% filter(SZ_PRS_AVG == 1)

Bottom_top_decile_merged <- merge(Bottom_decile_PRS_only, Top_decile_PES, by = "IID")

table(Bottom_top_decile_merged$SZ.y, Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)

SZ_dec <- glm(SZ.x ~ Sex.x + Age.x + Binary_decile, family = "binomial", data = Bottom_top_decile_merged)


