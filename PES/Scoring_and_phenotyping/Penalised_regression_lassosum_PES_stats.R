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
  make_option("--covar_file", action="store", default=NA, type='character',
              help="Covariate file"),
  make_option("--pheno_file", action="store", default=NA, type='character',
              help="Phenotype file"),
  make_option("--pathway", action="store", default=NA, type='character',
              help="Name of pathway"),
  make_option("--bfile_prefix", action="store", default=NA, type='character',
              help="Prefix for binary plink file-set")
)
opt = parse_args(OptionParser(option_list=option_list))

## Set system working directory

setwd("~/data/users/william/SZ_and_BIP_PES")

## Invoke 6 threads

c1 <- makeCluster(3)

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

Covariate_file <- Covariate_file %>% select(IID, Sex, Age, PC1, PC2, PC3, PC4, PC5, Batch)

cat("#########################")
cat("\n")
cat("Loading phenotype file")
cat("\n")
cat("#########################")
cat("\n")

Pheno_file <- fread(paste("",opt$pheno_file, sep=""), header = T)

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

cor <- p2cor(p = Sum_stats$P, n = 161405, sign=log(Sum_stats$OR))

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
                         A1=Sum_stats$A1, A2=Sum_stats$A2,
                         ref.bfile=ref.bfile, test.bfile=ref.bfile, 
                         LDblocks = LDblocks, cluster=c1)


cat("#########################")
cat("\n")
cat("lassosum finished, moving onto validation")
cat("\n")
cat("#########################")
cat("\n")


sv <- splitvalidate(out, pheno = Pheno_file, covar = Covariate_file)

file = paste(opt$pathway, '_lassosum', opt$phenotype, '.txt', sep="")

sink(file, append=T)

sv

sink()

#Clear environment after run
rm(list = ls())
