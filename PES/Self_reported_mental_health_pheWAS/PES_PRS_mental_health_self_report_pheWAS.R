#############################

## Correlation of PES/PRS with self-reported MHQ phenotypes

## SZ and BIP excluded + controls in the discovery analyses

## William Reay (2021)

############################

suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))

##Specify command line inputs
option_list = list(
  make_option("--disorder_name", action="store", default=NA, type='character',
              help="Name of the disorder for the PGS or PES [required]")
)
opt = parse_args(OptionParser(option_list=option_list))

print(opt)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_self_reported_mental_illness/")

## Read in PES/PRS

PES_PRS <- fread("../UKBB_MHQ_no_self_reported_biochem/Entire_EUR_genotyped_UKBB_scores/FINAL_SCORES_PES_PRS/All_SZ_BIP_PES_PRS_all_UKBB.txt", header = T)

PES_PRS <- PES_PRS[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28)]

## Load in MHQ self-report data

MHQ_self_report <- fread("UKBB_self_report_mental_illness_df_SZ_BIP_excl.txt", header = T)

MHQ_self_report <- rename(MHQ_self_report, "IID"="f.eid")

## Merge

Merged_MHQ_scoring <- merge(PES_PRS, MHQ_self_report, by = "IID")

## Define list of scores

Scores_to_test <- as.list(colnames(Merged_MHQ_scoring[,2:15]))

## Define MH col names

Codes <- fread("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_self_reported_mental_illness/Self_report_MH_coding.txt", header = T)

Name_input <- as.list(Codes$Input_meaning)

## Define function to automate regression

MHQ_self_report_mental_illness <- function(pheno_col, score_col, MH_col, df) {
  test_df <- df[df[[pheno_col]] == 1 | df[[MH_col]] == 0,]
  test_df$scaled_score <- as.numeric(scale(test_df[[score_col]]))
  fmla <- as.formula(paste(pheno_col, "~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 +  scaled_score + Batch"))
  mod <- glm(fmla, family = "binomial", data = test_df)
  return(summary(mod))
}

Run_MHQ <- sapply(Scores_to_test, MHQ_self_report_mental_illness, pheno_col =opt$disorder_name, MH_col = "No_mental_health", 
                  df = Merged_MHQ_scoring)

Run_MHQ_extract <- apply(Run_MHQ, 2, function(x) return(as.data.frame(x$coefficients)[15, 1:4]))

Output <- data.frame()
for (i in 1:length(Run_MHQ_extract)) {
  Output <- rbind(Output, Run_MHQ_extract[[i]])
}
rownames(Output) <- Scores_to_test

Output$FDR <- p.adjust(Output$'Pr(>|z|)', method="fdr")

Output$Disorder <- opt$disorder_name


write.table(Output, file = paste("Results/UKBB_", opt$disorder_name, "_results_PES_PGS", 
                                 ".txt", sep=""),
            sep = "\t", row.names = T, quote = F)

#Clear environment after run
rm(list = ls())
