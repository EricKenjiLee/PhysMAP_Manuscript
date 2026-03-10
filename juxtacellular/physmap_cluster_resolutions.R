###############################################################################
# physmap_cluster_resolutions.R
#
# Generates PhysMAP (WNN) UMAP embedding visualizations with cluster
# assignments color-coded at Leiden resolutions: 0.5, 1.0, 1.5, 2.0, 2.5.
#
# Uses k.nn = 20 for all neighbor computations.
#
# Output: combined PDF panel to juxtacellular/output/
###############################################################################

library(here)
library(ggExtra)
library(ggpubr)
library(scatterpie)
library(reticulate)
library(mclust)
library(patchwork)

set.seed(42)

here::i_am("README.md")

source(here::here("constants.R"))
source(here::here("juxtacellular", "helperFunctions.R"))

OUTPUT_DIR <- here::here("juxtacellular", "output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Parameters ---
pcs <- 1:30
nDims <- 30
numComponents <- 30
k.nn <- 20
resolutions <- c(0.5, 1.0, 1.5, 2.0, 2.5)

# --- Load data ---
allData <- readJianingData(here::here("juxtacellular", "JianingData", "MergedData.mat"))
juxtaData <- allData$data
juxtaData$isExcitatory <- ifelse(grepl("^E", trimws(juxtaData$CellType)),
                                 "Excitatory", "Inhibitory")

# --- Compute individual modality representations ---
# (required to produce the PCA reductions that FindMultiModalNeighbors needs)
tempFeat <- calcRepresentation(juxtaData, 'features', 5, 1:5, k.nn = k.nn)
juxtaData <- tempFeat$data

tempISI <- calcRepresentation(juxtaData, 'ISI', nDims, pcs, "ISI", TRUE,
                              nc = numComponents, metric = UMAP.metric,
                              k.nn = k.nn)
juxtaData <- tempISI$data

tempWF <- calcRepresentation(juxtaData, 'WF', nDims, pcs, "WF", FALSE,
                             nc = numComponents, metric = UMAP.metric,
                             k.nn = k.nn)
juxtaData <- tempWF$data

tempPSTH <- calcRepresentation(juxtaData, 'PSTH', nDims, pcs, "PSTH", TRUE,
                               nc = numComponents, metric = UMAP.metric,
                               k.nn = k.nn)
juxtaData <- tempPSTH$data

# --- Compute WNN (Pooled multimodal) embedding ---
juxtaData <- FindMultiModalNeighbors(
  juxtaData,
  reduction.list = list("WFpca", "ISIpca", "PSTHpca"),
  dims.list = list(1:30, 1:30, 1:30),
  k.nn = k.nn
)
juxtaData <- RunUMAP(juxtaData, nn.name = "weighted.nn",
                     reduction.name = "wnn.umap",
                     reduction.key = "wnnUMAP_",
                     seed.use = UMAP.SEED,
                     metric = "euclidean",
                     n.neighbors = 2)

# --- Generate cluster-colored plots at each resolution ---
plotList <- list()
clustersList <- list()

for (resV in resolutions) {
  juxtaData <- FindClusters(juxtaData, graph.name = "wsnn",
                            algorithm = ALGORITHM,
                            resolution = resV, verbose = FALSE)
  clustersList[[as.character(resV)]] <- FindClusters(juxtaData, graph.name = "wsnn",
                                                     algorithm = ALGORITHM,
                                                     resolution = resV, verbose = FALSE)

  wnnE <- as.data.frame(Embeddings(juxtaData, reduction = 'wnn.umap'))
  wnnE$cluster <- juxtaData$seurat_clusters
  wnnE$isExcitatory <- juxtaData$isExcitatory

  p <- ggplot(wnnE, aes(x = wnnUMAP_1, y = wnnUMAP_2,
                         color = cluster, shape = isExcitatory)) +
    geom_point(size = 2) +
    scale_shape_manual(values = c("Excitatory" = 16, "Inhibitory" = 16)) +
    theme_void() +
    ggtitle(paste("Resolution =", resV))

  plotList[[as.character(resV)]] <- p
}

# --- Combine into a single figure ---
pCombined <- plotList[[1]] | plotList[[2]] | plotList[[3]] | plotList[[4]] | plotList[[5]]

print(pCombined)
ggsave(file.path(OUTPUT_DIR, "PhysMAP_cluster_resolutions.pdf"),
       pCombined, width = 20, height = 4)
