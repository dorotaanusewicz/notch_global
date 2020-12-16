setwd("~/Desktop/monocle_new")

library(monocle3)
library(devtools)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)


express_matrix1 =  read.table("express_matrix1.txt", header = TRUE, 
                              sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE)
express_matrix2 =  read.table("express_matrix2.txt", header = TRUE, 
                              sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE)
express_matrix3 =  read.table("express_matrix3.txt", header = TRUE, 
                              sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE)
express_matrix4 =  read.table("express_matrix4.txt", header = TRUE, 
                              sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE)
express_matrix5 =  read.table("express_matrix5.txt", header = TRUE, 
                              sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE)
express_matrix6 =  read.table("express_matrix6.txt", header = TRUE, 
                              sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE)
express_matrix7 =  read.table("express_matrix7.txt", header = TRUE, 
                              sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE)


expression_matrix = cbind(express_matrix1, express_matrix2, express_matrix3,
                          express_matrix4, express_matrix5, express_matrix6, 
                          express_matrix7)

cell_metadata = read.table("cell_metadata.txt", header = TRUE, 
                           sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE)
gene_metadata = read.table("gene_metadata_NE.txt", header = TRUE, 
                           sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE)


rownames(expression_matrix) = expression_matrix$gene_short_name
row.names(cell_metadata) = colnames(expression_matrix)
row.names(gene_metadata) = row.names(expression_matrix)

expression_matrix = expression_matrix[,-1]
expression_matrix = as.matrix(expression_matrix)

row.names(cell_metadata) = colnames(expression_matrix)
row.names(gene_metadata) = row.names(expression_matrix)

library(Matrix)
exp_matric = as(expression_matrix, "dgCMatrix")

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)



cds2 <- preprocess_cds(cds2, num_dim = 100)
plot_pc_variance_explained(cds2) 

set.seed(4837)
cds2 <- reduce_dimension(cds2, reduction_method = "UMAP", preprocess_method = "PCA")

plot_cells(cds2)
plot_cells(cds2, color_cells_by="cancertype", cell_size = 1, group_label_size = 3.5)
plot_cells(cds, color_cells_by="cancersubtype")
plot_cells(cds, color_cells_by="tissuetype")