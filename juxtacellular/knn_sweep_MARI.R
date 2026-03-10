###############################################################################
# knn_sweep_MARI.R
#
# Sweeps k.nn = {5, 10, 20, 30, 40} and for each value:
#   1. Recomputes individual modality representations (ISI, WF, PSTH, Concat)
#   2. Recomputes multimodal WNN (Pooled) representation
#   3. Sweeps Leiden resolution 0.1--3.0, computing MARI at each resolution
#   4. Generates ARIg1 plot (MARI vs resolution, lines per modality)
#
# Outputs: PDF plots and CSV source data to juxtacellular/output/
###############################################################################

library(here)
library(ggExtra)
library(ggpubr)
library(scatterpie)
library(reticulate)
library(mclust)
library(reshape2)

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
knnValues <- c(5, 10, 20, 30, 50)
xV <- seq(0.1, 3.0, 0.1)

# --- Load data once ---
allData <- readJianingData(here::here("juxtacellular", "JianingData", "MergedData.mat"))
baseData <- allData$data
baseData$isExcitatory <- ifelse(grepl("^E", trimws(baseData$CellType)),
                                "Excitatory", "Inhibitory")

# --- Collector for combined results across all k.nn values ---
allResultsDF <- data.frame()

# --- Main k.nn sweep loop ---
for (knn in knnValues) {

  message(paste("=== Processing k.nn =", knn, "==="))

  # Start from a fresh copy of the base data for each k.nn
  juxtaData <- baseData

  # Step 1: Compute features representation
  tempFeat <- calcRepresentation(juxtaData, 'features', 5, 1:5, k.nn = knn)
  juxtaData <- tempFeat$data

  # Step 2: Compute ISI representation
  tempISI <- calcRepresentation(juxtaData, 'ISI', nDims, pcs, "ISI", TRUE,
                                nc = numComponents, metric = UMAP.metric,
                                k.nn = knn)
  juxtaData <- tempISI$data
  ARIvCell <- tempISI$layerCellTypeARI

  # Step 3: Compute WF representation
  tempWF <- calcRepresentation(juxtaData, 'WF', nDims, pcs, "WF", FALSE,
                               nc = numComponents, metric = UMAP.metric,
                               k.nn = knn)
  juxtaData <- tempWF$data
  ARIvCell <- cbind(ARIvCell, tempWF$layerCellTypeARI)

  # Step 4: Compute PSTH representation
  tempPSTH <- calcRepresentation(juxtaData, 'PSTH', nDims, pcs, "PSTH", TRUE,
                                 nc = numComponents, metric = UMAP.metric,
                                 k.nn = knn)
  juxtaData <- tempPSTH$data
  ARIvCell <- cbind(ARIvCell, tempPSTH$layerCellTypeARI)

  # Step 5: Compute Concat representation
  tempConcat <- calcRepresentation(juxtaData, 'concat', nDims, pcs, "concat", TRUE,
                                   nc = numComponents, metric = UMAP.metric,
                                   k.nn = knn)
  juxtaData <- tempConcat$data
  ARIvCell <- cbind(ARIvCell, tempConcat$layerCellTypeARI)

  # Step 6: Compute Multimodal WNN (Pooled)
  juxtaData <- FindMultiModalNeighbors(
    juxtaData,
    reduction.list = list("WFpca", "ISIpca", "PSTHpca"),
    dims.list = list(1:30, 1:30, 1:30),
    k.nn = knn
  )
  juxtaData <- RunUMAP(juxtaData, nn.name = "weighted.nn",
                       reduction.name = "wnn.umap",
                       reduction.key = "wnnUMAP_",
                       seed.use = UMAP.SEED,
                       metric = "euclidean",
                       n.neighbors = 2)

  # Sweep resolution for Pooled (WNN) representation
  d2_pooled <- c()
  for (resV in xV) {
    data <- FindClusters(juxtaData, graph.name = "wsnn",
                         algorithm = ALGORITHM,
                         resolution = resV, verbose = FALSE)
    d2_pooled <- c(d2_pooled, MARI(data$seurat_clusters, data$layerCellType))
  }

  # Step 7: Assemble ARIvCell and generate ARIg1 plot
  ARIvCell <- cbind(ARIvCell, d2_pooled, xV)
  ARIv <- data.frame(ARIvCell)
  colnames(ARIv) <- c('ISI', 'WF', 'PSTH', 'Concat', 'Pooled', 'xV')
  ARIvmelted <- melt(ARIv, id.vars = "xV")

  ARIg1 <- ggplot(data = ARIvmelted,
                  aes(x = xV, y = value, group = variable, col = variable)) +
    geom_line(linewidth = 2) +
    theme_minimal() +
    ggtitle(paste("MARI vs Resolution (k.nn =", knn, ")")) +
    xlab("Resolution") + ylab("MARI")

  print(ARIg1)
  ggsave(file.path(OUTPUT_DIR, paste0("ARIg1_knn_", knn, ".pdf")),
         ARIg1, width = 8, height = 6)

  # Accumulate results for combined analysis
  ARIvmelted$k.nn <- knn
  allResultsDF <- rbind(allResultsDF, ARIvmelted)
}

# --- Save combined CSV ---
write.csv(allResultsDF, file.path(OUTPUT_DIR, "MARI_knn_sweep_all.csv"),
          row.names = FALSE)
