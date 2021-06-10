##############################################

## Definining genes for input into the PES from each pathway

## William Reay (2021)

##############################################

##Load dependencies
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))

##Specify command line inputs
option_list = list(
  make_option("--gene_list_file", action="store", default=NA, type='character',
              help="The name of the gene-list to extract the genes [required]"),
  make_option("--pathway_name", action="store", default=NA, type='character',
              help="Name of the pathway")
)
opt = parse_args(OptionParser(option_list=option_list))

## Read in Autosomal genes with MHC removed

MHC_removed <- fread("~/cloudstor/Cross_disorder_PES_2019/MHC_removed_autosomes.txt", header = F)

MHC_removed <- rename(MHC_removed, "Gene"="V1", "CHR"="V2", "START"="V3","END"="V4","STRAND"="V5")

## Read in gene list

Gene_list <- read.table(file = paste("./",opt$gene_list_file,sep=""),
                                                 header=FALSE,sep="\t")

Gene_list <- rename(Gene_list, "Gene"="V1")

Merged <- merge(MHC_removed, Gene_list, by="Gene")

Merged <- Merged %>% select(CHR, START, END, STRAND)

write.table(Merged, file=paste(opt$pathway_name, '_genes_LIST.txt', sep = ""), sep="\t",
            row.names = F, quote = F)