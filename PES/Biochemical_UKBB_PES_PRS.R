################################

## Effect of PES/PRS on blood/urine biochemical traits amongst UKBB participants who did not self-report mental illness in the MHQ

## William Reay (2021)

#################################

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/")


##Specify command line inputs
option_list = list(
  make_option("--biochem_name", action="store", default=NA, type='character',
              help="The name of biochemical that will be analysed [required]"),
  make_option("--biochem_name_col_id", action="store", default=NA, type='character',
              help="UKBB column id for the biochemical trait [required]"),
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
cat("Constructing model to test the association of",  paste(opt$biochem_name, ' with ', opt$score_name, sep = ""))
cat("\n")
cat("#########################")
cat("\n")

## Read in df

Merged_biochem_data <- fread("")

cat("\n")
cat("#########################")
cat("\n")
cat("Dataframe imported with scores and biochemical data")
cat("\n")
cat("#########################")
cat("\n")

Biochem_PES_PRS <- function(biochem_col, score_col, df) {
  biochem_df <- df[!is.na(df$biochem_col),]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  mod <- lm(biochem_col ~ Sex*Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + score + Batch, data = seropos_df)
  return(summary(mod))
}


Run_lm <- Biochem_PES_PRS(opt$biochem_name_col_id, opt$score_col_id, Merged_biochem_data)

Output <- as.data.frame(Run_lm$coefficients)[16, 1:4]

write.table(Output, file = paste("UKBB_MHQ_no_self_reported_biochem/", opt$disorder_name, "_", opt$biochem_name, "_", 
                                 opt$score_name, ".txt", sep=""),
            sep = "\t", row.names = F, quote = F)

## Sex stratified


Biochem_PES_PRS_sex_strat <- function(biochem_col, score_col, df, sex) {
  biochem_df <- df[!is.na(df$biochem_col),]
  biochem_df <- biochem_df[biochem_df[[Sex]] == sex,]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  mod <- lm(biochem_col ~ Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + score + Batch, data = seropos_df)
  return(summary(mod))
}

Male_run_lm <- Biochem_PES_PRS(opt$biochem_name_col_id, opt$score_col_id, Merged_biochem_data, "Male")

Male_output <- as.data.frame(Male_run_lm$coefficients)[14, 1:4]

write.table(Male_output, file = paste("UKBB_MHQ_no_self_reported_biochem/Sex_stratified/Male_", opt$disorder_name, "_", opt$biochem_name, "_", 
                                 opt$score_name, ".txt", sep=""),
            sep = "\t", row.names = F, quote = F)

Female_run_lm <- Biochem_PES_PRS(opt$biochem_name_col_id, opt$score_col_id, Merged_biochem_data, "Female")

Female_output <- as.data.frame(Female_run_lm$coefficients)[14, 1:4]

write.table(Female_output, file = paste("UKBB_MHQ_no_self_reported_biochem/Sex_stratified/Female_", opt$disorder_name, "_", opt$biochem_name, "_", 
                                      opt$score_name, ".txt", sep=""),
            sep = "\t", row.names = F, quote = F)

#Clear environment after run
rm(list = ls())