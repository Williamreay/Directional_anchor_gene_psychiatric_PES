#############################################

## Association of PES and PRS with BIP in the UKBB cohort - VALIDATION

## BIP cases and double the number of controls (randomly selected from MHQ participants)

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

All_BIP_scores <- fread("Best_DA_gene_PES/BIP_best_con_or_lib/UKBB_BIP_VALIDATION/Combined_UKBB_BIP_replication_scores.txt", header = T, sep = "\t")

Colnames_all_BIP_scores <- make.unique(names(All_BIP_scores))

colnames(All_BIP_scores) <- Colnames_all_BIP_scores

## Read in phenotype and covariate data

BIP_pheno <- fread("UKBB_phenotypes/BIP_phenotypes/BIP_TEST_UKBB_pheno.txt", header = T)

BIP_cov <- fread("UKBB_phenotypes/BIP_phenotypes/Covariates_TEST_BIP_UKBB.txt", header = T)

## Merge all three df

BIP_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                     list(All_BIP_scores, BIP_pheno, BIP_cov))

BIP_merged$Batch <- as.factor(BIP_merged$Batch)

## Scale PES and PRS to have mean zero and unit variance (s.d. = 1)

BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)] <- lapply(BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)], function(x) c(scale(x)))

## Define function to test each score by itself and extract beta, SE, and P

## List of sets to test

Scores_to_test <- as.list(colnames(BIP_merged[,c(5, 10, 15, 20, 25, 30, 35)]))

## Function

PES_PRS_BIP_UKBB <- function(v){
  Score_model <- glm(glue::glue('BIP ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + {v} + Batch'), family = "binomial", data = BIP_merged)
  return(summary(Score_model)) 
}

## Iterate function over list of genes

BIP_score <- sapply(Scores_to_test, PES_PRS_BIP_UKBB)

## Extract beta, se, z, and p value for each gene

BIP_extract <- apply(BIP_score, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

BIP_results <- data.frame()

for (i in 1:length(BIP_extract)) {
  BIP_results <- rbind(BIP_results, BIP_extract[[i]])
}
rownames(BIP_results) <- Scores_to_test

write.csv(BIP_results, file="Best_DA_gene_PES/BIP_best_con_or_lib/UKBB_BIP_VALIDATION/Results/VALIDATION_BIP_UKBB_best_PES_PRS_results.csv", quote = F, row.names = T)


## Get variance explained on the liability scale (0.7% prevalence)

BIP_PES_PRS_r2 <- function(v, k, p){
  Null_model <- glm(BIP ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch, family = "binomial", data = BIP_merged)
  Score_model <- glm(glue::glue('BIP ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + Batch + {v}'), family = "binomial", data = BIP_merged)
  R2 <- nagelkerke(Score_model, null = Null_model)
  ## c.f. https://github.com/kn3in/ABC/blob/master/functions.R
  x <- qnorm(1 - k)
  z <- dnorm(x)
  i <- z / k
  cc <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
  theta <- i * ((p - k)/(1 - k)) * (i * ((p - k) / ( 1 - k)) - x)
  e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
  R2_liab <- cc * e * R2$Pseudo.R.squared.for.model.vs.null[3, ] / (1 + cc * e* theta * R2$Pseudo.R.squared.for.model.vs.null[3, ])
  return(R2_liab)
}

## Iterate function over list of genes

BIP_r2 <- sapply(Scores_to_test, BIP_PES_PRS_r2, k=0.01, p=0.33)

R2_liability_results <- data.frame()

for (i in 1:length(BIP_r2)) {
  R2_liability_results <- rbind(R2_liability_results, BIP_r2[[i]])
}
rownames(R2_liability_results) <- Scores_to_test

write.csv(R2_liability_results, file="Best_DA_gene_PES/BIP_best_con_or_lib/UKBB_BIP_VALIDATION/Results/VALIDATION_BIP_Liability_scale_r2_PES_PRS.csv", row.names = T, quote = F)

## Make forest plot of BIP OR per SD in the training and validation set

BIP_results$Score <- Scores_to_test

BIP_results$Cohort <- "Validation"

## Read in training results for best PES

BIP_train <- read_excel("Best_DA_gene_PES/BIP_best_PES_TRAIN_summary.xlsx")

BIP_train <- BIP_train %>% select(Score, Beta, SE)

BIP_train$Cohort <- "Training"

BIP_results <- BIP_results %>% select(Score, Estimate, 'Std. Error', Cohort)

BIP_results <- rename(BIP_results, "Beta"="Estimate", "SE"="Std. Error")

Plot_df <- rbind(BIP_train, BIP_results)


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

FP_BIP <- ggplot(data = Plot_df, aes(x=Score_name, y=OR, ymin=LOR, ymax=UOR, colour=Score_name)) +
  geom_pointrange() +
  geom_hline(yintercept=1, lty=2) +
  coord_flip() +
  ylab("Bipolar disorder OR [95% CI] per SD increase in PES") +
  theme_bw() +
  theme(legend.position = "null", axis.title.y = element_blank()) +
  facet_wrap(~Cohort,strip.position="top",nrow=2,scales = "free_y")

FP_BIP + scale_color_brewer(palette = "Dark2")

## Identify number of cases with top decile PES

Decile_df <- BIP_merged %>%
  mutate_at(vars(ends_with("AVG")), .funs = list(~ntile(.,10)))

## Assign individuals in the top PES decile

Top_decile_PES <- Decile_df %>%
  mutate_at(vars(ends_with("AVG")), .funs = list(~if_else(. == 10,1, 0)))

Top_decile_PES$PES_decile_sum <- rowSums(Top_decile_PES[, c(5, 10, 15, 20, 25, 35)])

Top_decile_PES$Binary_decile <- ifelse(Top_decile_PES$PES_decile_sum > 0, 1, 0)

table(Top_decile_PES$BIP, Top_decile_PES$PES_decile_sum)

Dec_tab <- table(Top_decile_PES$BIP, Top_decile_PES$Binary_decile)

Decile_BIP <- glm(BIP ~  Binary_decile + BIP_PRS_AVG + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
                    PC8 + PC9 + PC10 + Batch + Sex + Age,
                  family = "binomial", data = Top_decile_PES)

PRS_decile <- glm(Binary_decile ~ PRS,
                  family = "binomial", data = Top_decile_PES)

## Make mosaic plot of BIP cases and controls per decile grouping

mosaic( ~ BIP + Binary_decile, data = Top_decile_PES,
        highlighting = "BIP", highlighting_fill = c("grey", "plum"))  

Bottom_decile_PRS <- Decile_df %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~ntile(.,10))) %>%
  mutate_at(vars(starts_with("PRS")), .funs = list(~if_else(. == 1,1, 0)))

Bottom_decile_PRS_only <- Bottom_decile_PRS %>% filter(BIP_PRS_AVG == 1)

Bottom_top_decile_merged <- merge(Bottom_decile_PRS_only, Top_decile_PES, by = "IID")

table(Bottom_top_decile_merged$BIP.y, Bottom_top_decile_merged$PES_decile_sum)

Bottom_top_decile_merged$PES_decile_sum <- ifelse(Bottom_top_decile_merged$PES_decile_sum > 0, 1, 0)

BIP_dec <- glm(BIP.x ~ Sex.x + Age.x + Binary_decile, family = "binomial", data = Bottom_top_decile_merged)



