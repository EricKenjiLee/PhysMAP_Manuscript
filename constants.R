graphics.off()
library(Seurat)
library(cowplot)
library(dplyr)
library(R.matlab)
library(aricode)
library(ggplot2)

RESOLUTION = 2
ALGORITHM = 2
UMAP.SEED = 42

UMAP.neighbors = 20;
UMAP.mindist = 0.2
UMAP.metric = "cosine"
UMAP.components = 10
norm.margin = 4 # 1 is features, 2 is cells
