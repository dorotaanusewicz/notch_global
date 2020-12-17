setwd("~/Desktop/GLOBAL/heatmaps_pliki")

### receptor_ligand
eset = read.table("receptor_ligand_global_median.txt", sep = "\t", dec = ".", header = T)
rownames(eset) = eset$NAME
eset$NAME = NULL
eset = (as.matrix(eset))


hr <- hclust(as.dist(1-cor(t(eset), method="spearman")), method="complete")

library(gplots)

lighten <- function(color, factor = 3){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue = 255)
  col
}

col = greenred(50)
col = sapply(col, lighten)

# 1020 vs 654
heatmap.2(eset, Rowv=as.dendrogram(hr), Colv = NA, scale = "row", col = col, trace = "none", cexRow = 0.8, density.info = "none", dendrogram = "row",
          main = "Notch_ligands_receptors", xlab = "patients", ylab = "genes")


### modulators
mod = read.table("modulator_global_median.txt", sep = "\t", dec = ".", header = T)
rownames(mod) = mod$NAME
mod$NAME = NULL
mod = (as.matrix(mod))

hr <- hclust(as.dist(1-cor(t(mod), method="spearman")), method="complete")

library(gplots)

lighten <- function(color, factor = 3){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue = 255)
  col
}

col = greenred(50)
col = sapply(col, lighten)

# 1020 vs 654
heatmap.2(mod, Rowv=as.dendrogram(hr), Colv = NA, scale = "row", col = col, trace = "none", 
          density.info = "none", cexRow = 0.8, dendrogram = "row", main = "Notch_modulators", 
          xlab = "patients", ylab = "genes")

### transductors
trans = read.table("transductor_global_median.txt", sep = "\t", dec = ".", header = T)
rownames(trans) = trans$NAME
trans$NAME = NULL
trans = (as.matrix(trans))

hr <- hclust(as.dist(1-cor(t(trans), method="spearman")), method="complete")

library(gplots)

lighten <- function(color, factor = 3){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue = 255)
  col
}

col = greenred(50)
col = sapply(col, lighten)

# 1020 vs 654
heatmap.2(trans, Rowv=as.dendrogram(hr), Colv = NA, scale = "row", col = col, trace = "none", density.info = "none", 
          cexRow = 0.8, dendrogram = "row", main = "Notch_signal_transduction", 
          xlab = "patients", ylab = "genes")

### tft
tft = read.table("tft_global_median.txt", sep = "\t", dec = ".", header = T)
rownames(tft) = tft$NAME
tft$NAME = NULL
tft = (as.matrix(tft))

hr <- hclust(as.dist(1-cor(t(tft), method="spearman")), method="complete")

library(gplots)

lighten <- function(color, factor = 3){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue = 255)
  col
}

col = greenred(50)
col = sapply(col, lighten)

# 1020 vs 654
heatmap.2(tft, Rowv=as.dendrogram(hr), Colv = NA, scale = "row", col = col, trace = "none", density.info = "none", 
          cexRow = 0.8, main = "Notch_transcription factor", 
          xlab = "patients", ylab = "genes")
