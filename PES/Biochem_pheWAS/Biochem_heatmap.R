#####################

## Biochemical heatmap

## William Reay (2021)

######################

library(dplyr)
library(data.table)
library(ComplexHeatmap)
library(reshape2)
library(RColorBrewer)
library(readxl)
library(circlize)

setwd("~/Desktop/SZ_PES_mtCOJO_norm/Best_DA_gene_PES/UKBB_MHQ_no_self_reported_biochem/")

Hm_input <- read_excel("Combined_results_final_biochem_pheWAS.xlsx")

Biochem <- fread("Biochem_names.txt", header = T)

Biochem <- Biochem %>% filter(Biochem != "Oestradiol")

Hm_input <- merge(Hm_input, Biochem, by = "field_ID")

Hm_input <- Hm_input %>% select(t, Score, Biochem_trait.x)

Input <- dcast(Hm_input, Biochem_trait.x ~ Score, value.var = "t")

Input_mat <- as.matrix(Input)

rownames(Input_mat) <- Input_mat[, 1]

Input_mat <- Input_mat[, -1]

colnames(Input_mat) <- c("BIP CACNA1C PES", "BIP FADS1 PES",
                         "BIP FES PES", "BIP GRIN2A PES", "BIP PCCB PES",
                         "BIP PRS", "BIP RPS17 PES", "SZ CACNA1C PES",
                         "SZ FADS1 PES", "SZ FES PES", "SZ GRIN2A PES",
                         "SZ PCCB PES", "SZ PRS", "SZ RPS17 PES")


Input_mat <- `dimnames<-`(`dim<-`(as.numeric(Input_mat), dim(Input_mat)), dimnames(Input_mat))


col_fun = colorRamp2(c(-10, 0, 10), c("purple", "white", "green"))
col_fun(seq(-4, 4))

Heatmap(Input_mat, name = "t", rect_gp = gpar(col = "white", lwd = 2), 
       show_column_dend = TRUE, show_row_dend = TRUE, clustering_distance_rows = "pearson",
       row_dend_width = unit(1, "cm"), col = col_fun,
       row_names_gp = grid::gpar(fontsize = 8.5), column_names_gp = grid::gpar(fontsize = 7.5),
       column_names_rot = 50, column_title_gp = grid::gpar(fontsize = 11.5, fontface="bold"),
       cell_fun = function(j, i, x, y, w, h, fill) {
         if(Input_mat[i, j] > 3.85 || Input_mat[i, j] < -3.85 ) {
           grid.text("***", x, y)
         }
         if((Input_mat[i, j] > 2.89 && Input_mat[i, j] < 3.85) || (Input_mat[i, j] < -2.89 && Input_mat[i, j] > -3.85)) {
           grid.text("**", x, y)
         }
         if((Input_mat[i, j] > 1.96 && Input_mat[i, j] < 2.89) || (Input_mat[i, j] < -1.96 && Input_mat[i, j] > -2.89)) {
           grid.text("*", x, y)
         }
       }
)
