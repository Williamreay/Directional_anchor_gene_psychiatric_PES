################################

## Effect of PES/PRS on blood/urine biochemical traits amongst UKBB participants who did not self-report mental illness in the MHQ

## William Reay (2021)

#################################

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))

setwd("~/Desktop/SZ_PES_mtCOJO_norm/PES/")


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

cat("\n")
cat("#########################")
cat("\n")
cat("Read in list of biochemical traits to test")
cat("\n")
cat("#########################")
cat("\n")

Biochem_traits <- fread("~/Desktop/SZ_PES_mtCOJO_norm/PES/UKBB_MHQ_no_self_reported_biochem/List_of_biochem_field_IDs_both_sexes.txt",
                        header = F)

Oestradiol_excl <- Biochem_traits %>% filter(V1 != 'f.30790.0.0')

Female <- as.list(Biochem_traits)

Male_and_all <- as.list(Oestradiol_excl)

cat("\n")
cat("#########################")
cat("\n")
cat("Building regression models - male and females")
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


Run_lm <- sapply(Male_and_all, Biochem_PES_PRS,  score_col=opt$score_col_id, df=Merged_biochem_data)

Run_lm_extract <- apply(Run_lm, 2, function(x) return(as.data.frame(x$coefficients)[16, 1:4]))

Output <- data.frame()
for (i in 1:length(Run_lm_extract)) {
  Output <- rbind(Output, Run_lm_extract[[i]])
}
rownames(Output) <- Male_and_all

Output$FDR <- p.adjust(Output$'Pr(>|t|)', method="fdr")


write.table(Output, file = paste("UKBB_MHQ_no_self_reported_biochem/", opt$disorder_name, "_", 
                                 opt$score_name, ".txt", sep=""),
            sep = "\t", row.names = T, quote = F)

## Sex stratified


Biochem_PES_PRS_sex_strat <- function(biochem_col, score_col, df, sex) {
  biochem_df <- df[!is.na(df$biochem_col),]
  biochem_df <- biochem_df[biochem_df[[Sex]] == sex,]
  biochem_df$scaled_score <- as.numeric(scale(biochem_df[[score_col]]))
  mod <- lm(biochem_col ~ Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + score + Batch, data = seropos_df)
  return(summary(mod))
}

Male_run_lm <- sapply(Male_and_all, Biochem_PES_PRS,  score_col=opt$score_col_id, df=Merged_biochem_data, sex="Male")

Male_lm_extract <- apply(Male_run_lm, 2, function(x) return(as.data.frame(x$coefficients)[14, 1:4]))

Male_output <- data.frame()
for (i in 1:length(Male_run_lm_extract)) {
  Male_output <- rbind(Male_output, Male_run_lm_extract[[i]])
}
rownames(Male_output) <- Male_and_all

Male_output$FDR <- p.adjust(Male_output$'Pr(>|t|)', method="fdr")

write.table(Male_output, file = paste("UKBB_MHQ_no_self_reported_biochem/Sex_stratified/Male_", opt$disorder_name, "_", 
                                 opt$score_name, ".txt", sep=""),
            sep = "\t", row.names = T, quote = F)

Female_run_lm <- sapply(Female, Biochem_PES_PRS,  score_col=opt$score_col_id, df=Merged_biochem_data, sex="Female")

Female_lm_extract <- apply(Female_run_lm, 2, function(x) return(as.data.frame(x$coefficients)[14, 1:4]))

Female_output <- data.frame()
for (i in 1:length(Female_run_lm_extract)) {
  Female_output <- rbind(Female_output, Female_run_lm_extract[[i]])
}
rownames(Female_output) <- Female

Female_output$FDR <- p.adjust(Female_output$'Pr(>|t|)', method="fdr")


write.table(Female_output, file = paste("UKBB_MHQ_no_self_reported_biochem/Sex_stratified/Female_", opt$disorder_name, "_", 
                                      opt$score_name, ".txt", sep=""),
            sep = "\t", row.names = T, quote = F)

#Clear environment after run
rm(list = ls())