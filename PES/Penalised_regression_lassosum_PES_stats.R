######################################

## Penalised regression pipeline (lassosum)

## William Reay (2021)

#########################################

## Load dependencies

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(methods))
suppressMessages(library(magrittr))
suppressMessages(library(parallel))
suppressMessages(library(optparse))
suppressMessages(library(lassosum))

##Specify command line inputs
option_list = list(
  make_option("--phenotype_name", action="store", default=NA, type='character',
              help="The name of the trait to be analysed"),
  make_option("--sumstats_file", action="store", default=NA, type='character',
              help="File name for summary stats for trait to be analysed"),
  make_option("--sumstats_file", action="store", default=NA, type='character',
              help="File name for summary stats for trait to be analysed"),
  make_option("--covar_file", action="store", default=NA, type='character',
              help="File name for summary stats for trait to be analysed"),
  make_option("--bfile_prefix", action="store", default=NA, type='character',
              help="Prefix for binary plink file-set"),
  make_option("--sample_size", action="store", default=NA, type='character',
              help="Sample size for the GWAS to be profiled")
)
opt = parse_args(OptionParser(option_list=option_list))

## Set system working directory

setwd(system.file("~/ukbb2/SZ_and_BIP_UKBB_subsets", package="lassosum"))

## Invoke 6 threads

c1 <- makeCluster(6)

## Read in summary statistics

cat("#########################")
cat("\n")
cat("Loading summary statistics")
cat("\n")
cat("#########################")
cat("\n")

Sum_stats <- fread(paste("",opt$sumstats_file, sep=""), header = T)

cat("#########################")
cat("\n")
cat("Loading covariate file")
cat("\n")
cat("#########################")
cat("\n")

Covariate_file <- fread(paste("",opt$covar_file, sep=""), header = T)

## Specify prefix of the plink binary file-set which will be used to train/test the model via pseudovalidation

ref.bfile <- paste("",opt$bfile_prefix, sep="")

LDblocks <- "EUR.hg19"

## Calculate P-value wise correlations using inbuilt lasso sum function

cat("#########################")
cat("\n")
cat("Calculating correlation matrix")
cat("\n")
cat("#########################")
cat("\n")

cor <- p2cor(p = Sum_stats$P, n = paste("", opt$sample_size, sep=""), sign=log(Sum_stats$OR))

cat("#########################")
cat("\n")
cat("Correlation matrix finished")
cat("\n")
cat("#########################")
cat("\n")


cat("#########################")
cat("\n")
cat("Beginning lassosum penalised regression")
cat("\n")
cat("#########################")
cat("\n")



out <- lassosum.pipeline(cor=cor, chr=Sum_stats$CHR, pos=Sum_stats$BP, 
                         A1=Sum_stats$A1, Aw=Sum_stats$A2,
                         ref.bfile=ref.bfile, test.bfile=ref.bfile, 
                         LDblocks = LDblocks, cluster=cl)


cat("#########################")
cat("\n")
cat("lassosum finished, moving onto validation")
cat("\n")
cat("#########################")
cat("\n")



#Clear environment after run
rm(list = ls())
