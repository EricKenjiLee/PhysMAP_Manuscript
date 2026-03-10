###############################################################################
# MARI_resolution_sweep.R
#
# MARI vs Leiden resolution sweep for the CellExplorer dataset.
# Mirrors the approach in juxtacellular/MARI_GMM_sweep.R (Figure 2C).
#
# Modalities: WF, ISI, ACG, Features, PhysMAP (pooled WNN)
# Dataset: CellExplorer extracellular recordings
#
# Produces:
#   - Line plot of MARI vs Leiden resolution per modality
#   - PDF output + Excel source data
###############################################################################

library(tidyverse)
library(R.matlab)
library(Seurat)
library(aricode)
library(here)
library(reshape2)
library(openxlsx)

here::i_am("README.md")

############### Parameters ###############

UMAP.SEED <- 42
ALGORITHM <- 1
numpcs <- 40
dimV <- 1:40
RESOLUTION <- 1
xV <- seq(0.1, 3.0, 0.1)

OUTPUT_DIR <- here::here("CellExplorerv2", "output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
wb <- createWorkbook()

############### Load Data ###############

DATA_DIR <- here::here("CellExplorerv2", "Data")

featureData <- readMat(here::here("CellExplorerv2", "Data", "features_cellExp_new_Nov2023.mat"))
features <- featureData$features

acgData <- readMat(here::here("CellExplorerv2", "Data", "CellExplorer_ACG.mat"))

WFce <- readMat(here::here("CellExplorerv2", "Data", "finalWaveforms.mat"))
X_waveform <- WFce$X
cType <- WFce$CellTypeNames

## Build Seurat object with WF assay
dataSize <- dim(X_waveform)
cellIds <- seq(1, dataSize[1])
rownames(X_waveform) <- cellIds
colnames(X_waveform) <- seq(1, dataSize[2])
data <- CreateSeuratObject(counts = t(X_waveform), assay = "WF")
data@meta.data <- cbind(data@meta.data, cType)

## Add ISI assay
load(here::here("CellExplorerv2", "Data", "isi_cellExp.Rda"))
X_ISI <- t(isi)
dataSize <- dim(X_ISI)
rownames(X_ISI) <- cellIds
colnames(X_ISI) <- seq(1, dataSize[2])
ISI_assay <- CreateAssayObject(counts = t(X_ISI))
data[["ISI1"]] <- ISI_assay

## Add ACG assay
acg <- acgData$acgW
X_ACG <- t(acg)
dataSize <- dim(X_ACG)
rownames(X_ACG) <- cellIds
colnames(X_ACG) <- seq(1, dataSize[2])
ACG_assay <- CreateAssayObject(counts = t(X_ACG))
data[["ACG"]] <- ACG_assay

## Add features assay
X_features <- features
dataSize <- dim(X_features)
rownames(X_features) <- cellIds
colnames(X_features) <- seq(1, dataSize[2])
features_assay <- CreateAssayObject(counts = t(X_features))
data[["features"]] <- features_assay

############### Helper Functions ###############

runAnalysis <- function(data, whichAssay, RESOLUTION = 0.75,
                        metricV = "cosine", normalize = FALSE,
                        numPCS = numpcs, DimV = dimV, nc = 10) {
  DefaultAssay(data) <- whichAssay
  if (normalize) {
    data <- NormalizeData(data, normalization.method = "CLR", margin = 2)
  }
  data <- ScaleData(data)
  data <- FindVariableFeatures(data)
  data <- RunPCA(data, verbose = FALSE,
                 reduction.name = paste0(whichAssay, "pca"), npcs = numPCS)
  data <- RunPCA(data, verbose = FALSE, npcs = numPCS)
  data <- FindNeighbors(data, dims = DimV)
  data <- RunUMAP(data, dims = DimV,
                  reduction.name = paste0(whichAssay, "umap"),
                  metric = metricV, n.components = nc, seed.use = UMAP.SEED)
  data <- RunUMAP(data, dims = DimV,
                  reduction.name = paste0(whichAssay, "plotumap"),
                  metric = metricV)
  data <- FindClusters(data, algorithm = ALGORITHM,
                       resolution = RESOLUTION, verbose = FALSE)
  return(data)
}

calcARI <- function(data, whichAssay, res = xV) {
  d <- c()
  DefaultAssay(data) <- whichAssay
  for (resV in res) {
    data <- FindClusters(data, algorithm = ALGORITHM,
                         resolution = resV, verbose = FALSE)
    d <- c(d, MARI(data$seurat_clusters, data$cType))
  }
  return(d)
}

############### Run Analysis Per Modality ###############

# Features (10 PCs)
data <- runAnalysis(data, "features", RESOLUTION = RESOLUTION,
                    normalize = TRUE, numPCS = 10, DimV = 1:10)
dF <- calcARI(data, "features")

# ACG (40 PCs)
data <- runAnalysis(data, "ACG", RESOLUTION = RESOLUTION, normalize = TRUE)
dA <- calcARI(data, "ACG")

# ISI (40 PCs)
data <- runAnalysis(data, "ISI1", RESOLUTION = RESOLUTION, normalize = TRUE)
dISI <- calcARI(data, "ISI1")

# WF (40 PCs)
data <- runAnalysis(data, "WF", RESOLUTION = RESOLUTION, normalize = TRUE)
dWF <- calcARI(data, "WF")

############### PhysMAP (Pooled WNN) ###############

data <- FindMultiModalNeighbors(
  data,
  reduction.list = list("WFpca", "ISI1pca", "ACGpca", "featurespca"),
  dims.list = list(1:40, 1:40, 1:40, 1:10)
)
data <- RunUMAP(data, nn.name = "weighted.nn",
                reduction.name = "wnn.umap",
                reduction.key = "wnnUMAP_", seed.use = UMAP.SEED)

dWNN <- c()
for (resV in xV) {
  data <- FindClusters(data, graph.name = "wsnn", algorithm = ALGORITHM,
                       resolution = resV, verbose = FALSE)
  dWNN <- c(dWNN, MARI(data$seurat_clusters, data$cType))
}

############### Build MARI DataFrame & Plot ###############

ariDF <- data.frame(
  resolution = rep(xV, 5),
  MARI = c(dWF, dISI, dA, dF, dWNN),
  Modality = rep(c("WF", "ISI", "ACG", "Features", "PhysMAP"),
                 each = length(xV))
)

figMARI <- ggplot(ariDF, aes(x = resolution, y = MARI, color = Modality)) +
  geom_line(linewidth = 1.5) +
  theme_minimal() +
  ggtitle("MARI vs Leiden Resolution (CellExplorer)") +
  xlab("Resolution") + ylab("MARI")
show(figMARI)

ggsave(file.path(OUTPUT_DIR, "CellExplorer_MARI_resolution.pdf"),
       figMARI, width = 8, height = 6)

## Save source data
addWorksheet(wb, "MARI_sweep")
writeData(wb, "MARI_sweep", ariDF)
saveWorkbook(wb, file.path(OUTPUT_DIR, "CellExplorer_MARI_source_data.xlsx"),
             overwrite = TRUE)
