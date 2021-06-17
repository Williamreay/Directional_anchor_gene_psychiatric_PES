################################

## Effect of PES/PRS on blood/urine biochemical traits amongst UKBB participants who did not self-report mental illness in the MHQ

## IRNT transformation

## William Reay (2021)

#################################

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(RNOmni))

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem/")


##Specify command line inputs
option_list = list(
  make_option("--score_name", action="store", default=NA, type='character',
              help="Name of PGS or PES [required]"),
  make_option("--score_col_id", action="store", default=NA, type='character',
              help="Column id for the PGS or PES [required]"),
  make_option("--disorder_name", action="store", default=NA, type='character',
              help="Name of the disorder for the PGS or PES [required]")
)
opt = parse_args(OptionParser(option_list=option_list))

print(opt)

cat("\n")
cat("#########################")
cat("\n")
cat("Constructing model to test the association with ", paste(opt$score_name), sep = "")
cat("\n")
cat("#########################")
cat("\n")

## Read in df

Merged_biochem_data <- fread("Merged_all_scores_all_biochem.txt")

Merged_biochem_data$Batch <- as.factor(Merged_biochem_data$Batch)

cat("\n")
cat("#########################")
cat("\n")
cat("Dataframe imported with scores and biochemical data")
cat("\n")
cat("#########################")
cat("\n")

cat("\n")
cat("#########################")
cat("\n")
cat("Read in list of biochemical traits to test")
cat("\n")
cat("#########################")
cat("\n")

Biochem_traits <- fread("List_of_biochem_field_IDs_both_sexes.txt",
                        header = F)

Oestradiol_excl <- Biochem_traits %>% filter(V1 != 'f.30790.0.0')

Female <- as.list(Biochem_traits$V1)

Male_and_all <- as.list(Oestradiol_excl$V1)

cat("\n")
cat("#########################")
cat("\n")
cat("Building regression models - male and females")
cat("\n")
cat("#########################")
cat("\n")

IRNT <- function(biochem_col, score_col, df) {
  biochem_df <- df[!is.na(df[[biochem_col]]),]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  fmla <- as.formula(paste(biochem_col," ~ Sex*Age + Age2"))
  mod <- lm(fmla, data = biochem_df)
  Biochem_resid <- as.numeric(mod$residuals)
  ## Blom transformation
  Biochem_residuals_NORMALISED <- as.data.frame(RankNorm(Biochem_resid), k=3/8)
  biochem_df$Norm_resid <- Biochem_residuals_NORMALISED$'RankNorm(Biochem_resid)'
  Resid_mod <- lm(Norm_resid ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
                    PC8 + PC9 + PC10 + scaled_score + Batch, data = biochem_df)
  return(summary(Resid_mod))
}


Run_lm <- sapply(Male_and_all, IRNT,  score_col=opt$score_col_id, df=Merged_biochem_data)

Run_lm_extract <- apply(Run_lm, 2, function(x) return(as.data.frame(x$coefficients)[12, 1:4]))

Output <- data.frame()
for (i in 1:length(Run_lm_extract)) {
  Output <- rbind(Output, Run_lm_extract[[i]])
}
rownames(Output) <- Male_and_all

Output$FDR <- p.adjust(Output$'Pr(>|t|)', method="fdr")
Output$Score <- opt$score_name
Output$Biochem <- unlist(Male_and_all)

write.table(Output, file = paste("IRNT_Combined_results/", opt$disorder_name, "_IRNT_", 
                                 opt$score_name, ".txt", sep=""),
            sep = "\t", row.names = T, quote = F)
