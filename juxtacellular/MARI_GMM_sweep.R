###############################################################################
# MARI_GMM_sweep.R
#
# Produces: Figures 1, 2, S2, S3, S4, S5, S6, S7
# Dataset:  Juxtacellular mouse S1 (Yu et al., n=246 neurons)
#
# Adapted to run from the PhysMAP_Manuscript repo structure.
###############################################################################

library(here)
here::i_am("README.md")

source(here::here("constants.R"))
source(here::here("juxtacellular", "helperFunctions.r"))
library(mclust)
library(ggpubr)
library(caret)
library(nnet)
library(stringr)
library(reshape2)
library(openxlsx)

OUTPUT_DIR <- here::here("juxtacellular", "output")
SOURCE_DATA_DIR <- here::here("juxtacellular", "output", "source_data")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
wb <- createWorkbook()

JUXTA_DATA <- here::here("juxtacellular", "JianingData")
pcs <- 1:30
nDims <- 30
UMAP.components <- 10  # Override: 10D UMAP for classification (Desktop version used 10)
numComponents <- UMAP.components

allData <- readJianingData(file.path(JUXTA_DATA, "MergedData.mat"))
juxtaData <- allData$data

############### Figure 1A: WF UMAP ###############

wfResult <- calcRepresentation(juxtaData, 'WF', nDims, pcs, "WF", FALSE,
                               nc = numComponents, metric = UMAP.metric)
juxtaData <- wfResult$data

fig1A <- DimPlot(juxtaData, reduction = 'umap', group.by = "layerCellType", pt.size = 2) +
  ggtitle("Waveform Shape") + theme_minimal()
print(fig1A)
ggsave(file.path(OUTPUT_DIR, "Figure_1A_WF_UMAP.pdf"), fig1A, width = 7, height = 7)
addWorksheet(wb, "Fig1A")
fig1A_data <- data.frame(Embeddings(juxtaData, "umap"), cellType = juxtaData$layerCellType)
colnames(fig1A_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "Fig1A", fig1A_data)

############### Figure 1B: ISI UMAP ###############

isiResult <- calcRepresentation(juxtaData, 'ISI', nDims, pcs, "ISI", TRUE,
                                nc = numComponents, metric = UMAP.metric)
juxtaData <- isiResult$data

fig1B <- DimPlot(juxtaData, reduction = 'umap', group.by = "layerCellType", pt.size = 2) +
  ggtitle("ISI Distribution") + theme_minimal()
print(fig1B)
ggsave(file.path(OUTPUT_DIR, "Figure_1B_ISI_UMAP.pdf"), fig1B, width = 7, height = 7)
addWorksheet(wb, "Fig1B")
fig1B_data <- data.frame(Embeddings(juxtaData, "umap"), cellType = juxtaData$layerCellType)
colnames(fig1B_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "Fig1B", fig1B_data)

############### Figure 1C: WNN Weights (2 modalities) ###############

juxtaData <- FindMultiModalNeighbors(
  juxtaData,
  reduction.list = list("WFpca", "ISIpca"),
  dims.list = list(pcs, pcs),
  modality.weight.name = c("WF.weight", "ISI.weight"),
  verbose = FALSE
)

wfWeights <- juxtaData$WF.weight
isiWeights <- juxtaData$ISI.weight

weightDF <- data.frame(weight = wfWeights, modality = "WF")
weightDF <- rbind(weightDF, data.frame(weight = isiWeights, modality = "ISI"))

fig1C_hist <- ggplot(weightDF, aes(x = weight, fill = modality)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() + ggtitle("WNN Modality Weights")
print(fig1C_hist)
ggsave(file.path(OUTPUT_DIR, "Figure_1C_WNN_Weights_Hist.pdf"), fig1C_hist, width = 7, height = 5)
addWorksheet(wb, "Fig1C_Hist")
writeData(wb, "Fig1C_Hist", weightDF)

overallWF <- mean(wfWeights)
overallISI <- mean(isiWeights)
pieDF <- data.frame(modality = c("WF", "ISI"),
                    proportion = c(overallWF, overallISI))
fig1C_pie <- ggplot(pieDF, aes(x = "", y = proportion, fill = modality)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") + theme_void() +
  ggtitle(paste0("WF: ", round(overallWF * 100), "% | ISI: ", round(overallISI * 100), "%"))
print(fig1C_pie)
ggsave(file.path(OUTPUT_DIR, "Figure_1C_WNN_Weights_Pie.pdf"), fig1C_pie, width = 5, height = 5)
addWorksheet(wb, "Fig1C_Pie")
writeData(wb, "Fig1C_Pie", pieDF)

umapCoords2d <- Embeddings(juxtaData, reduction = "WFumap2d")
scatterDF <- data.frame(UMAP1 = umapCoords2d[, 1], UMAP2 = umapCoords2d[, 2],
                        WF_weight = wfWeights)
fig1C_scatter <- ggplot(scatterDF, aes(x = UMAP1, y = UMAP2, color = WF_weight)) +
  geom_point(size = 2) + scale_color_viridis_c() +
  theme_minimal() + ggtitle("Per-unit WNN Weight")
print(fig1C_scatter)
ggsave(file.path(OUTPUT_DIR, "Figure_1C_WNN_Weights_Scatter.pdf"), fig1C_scatter, width = 7, height = 7)
addWorksheet(wb, "Fig1C_Scatter")
writeData(wb, "Fig1C_Scatter", scatterDF)

############### Figure 1D: PhysMAP WNN UMAP (2 modalities) ###############

juxtaData <- RunUMAP(juxtaData, nn.name = "weighted.nn",
                     reduction.name = "wnn.umap",
                     reduction.key = "wnnUMAP_",
                     seed.use = UMAP.SEED)

fig1D <- DimPlot(juxtaData, reduction = 'wnn.umap', group.by = "layerCellType", pt.size = 2) +
  ggtitle("PhysMAP (WNN)") + theme_minimal()
print(fig1D)
ggsave(file.path(OUTPUT_DIR, "Figure_1D_PhysMAP_WNN_UMAP.pdf"), fig1D, width = 7, height = 7)
addWorksheet(wb, "Fig1D")
fig1D_data <- data.frame(Embeddings(juxtaData, "wnn.umap"), cellType = juxtaData$layerCellType)
colnames(fig1D_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "Fig1D", fig1D_data)

############### Figure 1E: Concatenated UMAP ###############

concatResult <- calcRepresentation(juxtaData, 'concat', nDims, pcs, "Concatenated",
                                   FALSE, nc = numComponents, metric = UMAP.metric)
juxtaData <- concatResult$data

fig1E <- DimPlot(juxtaData, reduction = 'umap', group.by = "layerCellType", pt.size = 2) +
  ggtitle("Concatenated (Unweighted)") + theme_minimal()
print(fig1E)
ggsave(file.path(OUTPUT_DIR, "Figure_1E_Concat_UMAP.pdf"), fig1E, width = 7, height = 7)
addWorksheet(wb, "Fig1E")
fig1E_data <- data.frame(Embeddings(juxtaData, "umap"), cellType = juxtaData$layerCellType)
colnames(fig1E_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "Fig1E", fig1E_data)

############### Figure 2A: PhysMAP with marker size (spike width, P2T, latency) ###############

load(file.path(JUXTA_DATA, "width.Rda"))
load(file.path(JUXTA_DATA, "ratio_p2t.Rda"))

wnnEmbed <- Embeddings(juxtaData, reduction = "wnn.umap")
sizeDF <- data.frame(UMAP1 = wnnEmbed[, 1], UMAP2 = wnnEmbed[, 2],
                     width = as.numeric(juxtaData$width),
                     ratio_p2t = as.numeric(juxtaData$ratio_p2t),
                     latency = as.numeric(juxtaData$latency))

fig2A_width <- ggplot(sizeDF, aes(x = UMAP1, y = UMAP2, size = width)) +
  geom_point(alpha = 0.7) + scale_size_continuous(range = c(0.5, 5)) +
  theme_minimal() + ggtitle("Spike Width")
print(fig2A_width)
ggsave(file.path(OUTPUT_DIR, "Figure_2A_SpikeWidth.pdf"), fig2A_width, width = 7, height = 7)

fig2A_p2t <- ggplot(sizeDF, aes(x = UMAP1, y = UMAP2, size = ratio_p2t)) +
  geom_point(alpha = 0.7) + scale_size_continuous(range = c(0.5, 5)) +
  theme_minimal() + ggtitle("Peak-to-Trough Ratio")
print(fig2A_p2t)
ggsave(file.path(OUTPUT_DIR, "Figure_2A_P2T.pdf"), fig2A_p2t, width = 7, height = 7)

sizeDF$logLatency <- log(sizeDF$latency + 1)
fig2A_latency <- ggplot(sizeDF, aes(x = UMAP1, y = UMAP2, size = logLatency)) +
  geom_point(alpha = 0.7) + scale_size_continuous(range = c(0.5, 5)) +
  theme_minimal() + ggtitle("Onset Latency (log)")
print(fig2A_latency)
ggsave(file.path(OUTPUT_DIR, "Figure_2A_Latency.pdf"), fig2A_latency, width = 7, height = 7)
addWorksheet(wb, "Fig2A")
writeData(wb, "Fig2A", sizeDF)

############### Figure 2B: PhysMAP colored by cell type and Leiden clusters ###############

juxtaData <- FindNeighbors(juxtaData, reduction = "wnn.umap", dims = 1:2)
juxtaData <- FindClusters(juxtaData, algorithm = ALGORITHM, resolution = RESOLUTION, verbose = FALSE)

fig2B_celltype <- DimPlot(juxtaData, reduction = 'wnn.umap', group.by = "layerCellType", pt.size = 2) +
  ggtitle("Cell Type") + theme_minimal()
print(fig2B_celltype)
ggsave(file.path(OUTPUT_DIR, "Figure_2B_CellType.pdf"), fig2B_celltype, width = 7, height = 7)
addWorksheet(wb, "Fig2B_CellType")
fig2B_ct <- data.frame(Embeddings(juxtaData, "wnn.umap"), cellType = juxtaData$layerCellType)
colnames(fig2B_ct)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "Fig2B_CellType", fig2B_ct)

fig2B_clusters <- DimPlot(juxtaData, reduction = 'wnn.umap', group.by = "seurat_clusters", pt.size = 2) +
  ggtitle("Leiden Clusters (res=2)") + theme_minimal()
print(fig2B_clusters)
ggsave(file.path(OUTPUT_DIR, "Figure_2B_Clusters.pdf"), fig2B_clusters, width = 7, height = 7)
addWorksheet(wb, "Fig2B_Clusters")
fig2B_cl <- data.frame(Embeddings(juxtaData, "wnn.umap"), cluster = juxtaData$seurat_clusters)
colnames(fig2B_cl)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "Fig2B_Clusters", fig2B_cl)

############### Figure 2C: MARI scores vs Leiden resolution ###############

xV <- seq(0.1, 3.0, 0.1)

isiARI <- isiResult$cellTypeARI
wfARI <- wfResult$cellTypeARI
concatARI <- concatResult$cellTypeARI

pooledARI <- c()
for (resV in xV) {
  juxtaData <- FindClusters(juxtaData, algorithm = 3, resolution = resV, verbose = FALSE)
  pooledARI <- c(pooledARI, MARI(juxtaData$seurat_clusters, juxtaData$CellType))
}

psthResult <- calcRepresentation(juxtaData, 'PSTH', nDims, pcs, "PSTH", TRUE,
                                 nc = numComponents, metric = UMAP.metric)
juxtaData <- psthResult$data
psthARI <- psthResult$cellTypeARI

ariDF <- data.frame(
  resolution = rep(xV, 5),
  MARI = c(isiARI, wfARI, psthARI, concatARI, pooledARI),
  Modality = rep(c("ISI", "WF", "PSTH", "Concat", "PhysMAP"), each = length(xV))
)

fig2C <- ggplot(ariDF, aes(x = resolution, y = MARI, color = Modality)) +
  geom_line(linewidth = 1) +
  theme_minimal() + ggtitle("MARI vs Leiden Resolution") +
  xlab("Resolution") + ylab("MARI")
print(fig2C)
ggsave(file.path(OUTPUT_DIR, "Figure_2C_MARI.pdf"), fig2C, width = 8, height = 6)
addWorksheet(wb, "Fig2C")
writeData(wb, "Fig2C", ariDF)

############### Figure 2D: Pie charts per Leiden cluster ###############

juxtaData <- FindClusters(juxtaData, algorithm = ALGORITHM, resolution = RESOLUTION, verbose = FALSE)

clusterIDs <- sort(unique(juxtaData$seurat_clusters))
fig2D_collector <- list()
for (cl in clusterIDs) {
  idx <- which(juxtaData$seurat_clusters == cl)
  cellTypes <- juxtaData$layerCellType[idx]
  ctTable <- table(cellTypes) / length(cellTypes)
  ctDF <- data.frame(cellType = names(ctTable), proportion = as.numeric(ctTable))
  ctDF$cluster <- cl
  fig2D_collector[[length(fig2D_collector) + 1]] <- ctDF

  innerPie <- ggplot(ctDF, aes(x = "", y = proportion, fill = cellType)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") + theme_void() + ggtitle(paste("Cluster", cl, "- Cell Types"))
  print(innerPie)

  wfW <- mean(juxtaData$WF.weight[idx])
  isiW <- mean(juxtaData$ISI.weight[idx])
  modDF <- data.frame(modality = c("WF", "ISI"), weight = c(wfW, isiW))

  outerPie <- ggplot(modDF, aes(x = "", y = weight, fill = modality)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") + theme_void() + ggtitle(paste("Cluster", cl, "- Modality Weights"))
  print(outerPie)
}
addWorksheet(wb, "Fig2D")
writeData(wb, "Fig2D", do.call(rbind, fig2D_collector))

############### Figure 2E: PSTHs per Leiden cluster ###############

fig2E_collector <- list()
DefaultAssay(juxtaData) <- "PSTH"
psthMat <- as.matrix(GetAssayData(juxtaData, layer = "counts"))
psthMat <- t(psthMat)

for (cl in clusterIDs) {
  idx <- which(juxtaData$seurat_clusters == cl)
  clusterPSTH <- psthMat[idx, ]
  meanPSTH <- colMeans(clusterPSTH)
  smoothPSTH <- stats::filter(meanPSTH, rep(1 / 5, 5), sides = 2)
  smoothPSTH[is.na(smoothPSTH)] <- meanPSTH[is.na(smoothPSTH)]
  normPSTH <- (smoothPSTH - min(smoothPSTH)) / (max(smoothPSTH) - min(smoothPSTH) + 1e-10)

  psthPlotDF <- data.frame(time = seq_along(normPSTH), response = as.numeric(normPSTH))
  psthPlotDF$cluster <- cl
  fig2E_collector[[length(fig2E_collector) + 1]] <- psthPlotDF
  p <- ggplot(psthPlotDF, aes(x = time, y = response)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
    theme_minimal() + ggtitle(paste("Cluster", cl, "PSTH")) +
    xlab("Time (bins)") + ylab("Normalized Response")
  print(p)
}
addWorksheet(wb, "Fig2E")
writeData(wb, "Fig2E", do.call(rbind, fig2E_collector))

############### Figure 2F: Classification accuracy by cell type and modality ###############

# Compute high-dimensional (10D) WNN UMAP for classification
juxtaData <- RunUMAP(juxtaData, nn.name = "weighted.nn",
                     reduction.name = "wnn.umap.hd",
                     reduction.key = "wnnUMAPHD_",
                     seed.use = UMAP.SEED,
                     n.components = UMAP.components)

wnnEmbed2F <- data.frame(Embeddings(juxtaData, reduction = "wnn.umap.hd"))
wfEmbed2F <- data.frame(Embeddings(juxtaData, reduction = "WFumap"))
isiEmbed2F <- data.frame(Embeddings(juxtaData, reduction = "ISIumap"))
psthEmbed2F <- data.frame(Embeddings(juxtaData, reduction = "PSTHumap"))

featDF2F <- data.frame(width = as.numeric(juxtaData$width),
                       ratio_p2t = as.numeric(juxtaData$ratio_p2t))

allResults <- list()
modNames <- c("PhysMAP", "WF", "ISI", "PSTH", "Features")
embedList <- list(wnnEmbed2F, wfEmbed2F, isiEmbed2F, psthEmbed2F, featDF2F)

for (i in seq_along(modNames)) {
  res <- doClassifyJuxta(embedList[[i]], juxtaData, numreps = 5,
                         method = 'boot', seedV = 1, whichType = 'layercells')
  res$uF$Modality <- modNames[i]
  allResults[[i]] <- res$uF
}

classDF <- do.call(rbind, allResults)
classDF$cellClass <- gsub("Class: ", "", classDF$cellClass)

fig2F <- ggplot(classDF, aes(x = cellClass, y = AccV, color = Modality, group = Modality)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  theme_minimal() + ggtitle("Classification Accuracy by Cell Type") +
  xlab("Cell Type") + ylab("Balanced Accuracy (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(fig2F)
ggsave(file.path(OUTPUT_DIR, "Figure_2F_Classification.pdf"), fig2F, width = 9, height = 6)
addWorksheet(wb, "Fig2F")
writeData(wb, "Fig2F", classDF)

############### Figure S2A: WF UMAP (3-modality context) ###############

figS2A <- DimPlot(juxtaData, reduction = 'WFumap2d', group.by = "layerCellType", pt.size = 2) +
  ggtitle("WF UMAP") + theme_minimal()
print(figS2A)
ggsave(file.path(OUTPUT_DIR, "Figure_S2A_WF_UMAP.pdf"), figS2A, width = 7, height = 7)
addWorksheet(wb, "FigS2A")
figS2A_data <- data.frame(Embeddings(juxtaData, "WFumap2d"), cellType = juxtaData$layerCellType)
colnames(figS2A_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "FigS2A", figS2A_data)

############### Figure S2B: PSTH UMAP ###############

figS2B <- DimPlot(juxtaData, reduction = 'PSTHumap2d', group.by = "layerCellType", pt.size = 2) +
  ggtitle("PSTH UMAP") + theme_minimal()
print(figS2B)
ggsave(file.path(OUTPUT_DIR, "Figure_S2B_PSTH_UMAP.pdf"), figS2B, width = 7, height = 7)
addWorksheet(wb, "FigS2B")
figS2B_data <- data.frame(Embeddings(juxtaData, "PSTHumap2d"), cellType = juxtaData$layerCellType)
colnames(figS2B_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "FigS2B", figS2B_data)

############### Figure S2C: ISI UMAP (3-modality context) ###############

figS2C <- DimPlot(juxtaData, reduction = 'ISIumap2d', group.by = "layerCellType", pt.size = 2) +
  ggtitle("ISI UMAP") + theme_minimal()
print(figS2C)
ggsave(file.path(OUTPUT_DIR, "Figure_S2C_ISI_UMAP.pdf"), figS2C, width = 7, height = 7)
addWorksheet(wb, "FigS2C")
figS2C_data <- data.frame(Embeddings(juxtaData, "ISIumap2d"), cellType = juxtaData$layerCellType)
colnames(figS2C_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "FigS2C", figS2C_data)

############### Figure S2D: 3-modality WNN weights ###############

juxtaData <- FindMultiModalNeighbors(
  juxtaData,
  reduction.list = list("WFpca", "ISIpca", "PSTHpca"),
  dims.list = list(pcs, pcs, pcs),
  modality.weight.name = c("WF.weight.3", "ISI.weight.3", "PSTH.weight.3"),
  verbose = FALSE
)

wfW3 <- juxtaData$WF.weight.3
isiW3 <- juxtaData$ISI.weight.3
psthW3 <- juxtaData$PSTH.weight.3

weightDF3 <- data.frame(
  weight = c(wfW3, isiW3, psthW3),
  modality = rep(c("WF", "ISI", "PSTH"), each = ncol(juxtaData))
)

figS2D_hist <- ggplot(weightDF3, aes(x = weight, fill = modality)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() + ggtitle("3-Modality WNN Weights")
print(figS2D_hist)
ggsave(file.path(OUTPUT_DIR, "Figure_S2D_3Mod_Weights_Hist.pdf"), figS2D_hist, width = 7, height = 5)
addWorksheet(wb, "FigS2D_Hist")
writeData(wb, "FigS2D_Hist", weightDF3)

pieDF3 <- data.frame(
  modality = c("WF", "ISI", "PSTH"),
  proportion = c(mean(wfW3), mean(isiW3), mean(psthW3))
)
figS2D_pie <- ggplot(pieDF3, aes(x = "", y = proportion, fill = modality)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") + theme_void() +
  ggtitle(paste0("WF: ", round(mean(wfW3) * 100), "% | ISI: ",
                 round(mean(isiW3) * 100), "% | PSTH: ", round(mean(psthW3) * 100), "%"))
print(figS2D_pie)
ggsave(file.path(OUTPUT_DIR, "Figure_S2D_3Mod_Weights_Pie.pdf"), figS2D_pie, width = 5, height = 5)
addWorksheet(wb, "FigS2D_Pie")
writeData(wb, "FigS2D_Pie", pieDF3)

############### Figure S2E: 3-modality PhysMAP UMAP ###############

juxtaData <- RunUMAP(juxtaData, nn.name = "weighted.nn",
                     reduction.name = "wnn.umap.3mod",
                     reduction.key = "wnnUMAP3_",
                     seed.use = UMAP.SEED)

figS2E <- DimPlot(juxtaData, reduction = 'wnn.umap.3mod', group.by = "layerCellType", pt.size = 2) +
  ggtitle("PhysMAP (3 Modalities)") + theme_minimal()
print(figS2E)
ggsave(file.path(OUTPUT_DIR, "Figure_S2E_PhysMAP_3Mod.pdf"), figS2E, width = 7, height = 7)
addWorksheet(wb, "FigS2E")
figS2E_data <- data.frame(Embeddings(juxtaData, "wnn.umap.3mod"), cellType = juxtaData$layerCellType)
colnames(figS2E_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "FigS2E", figS2E_data)

############### Figure S2F: 3-modality concatenated UMAP ###############

figS2F <- DimPlot(juxtaData, reduction = 'concatumap2d', group.by = "layerCellType", pt.size = 2) +
  ggtitle("Concatenated (3 Modalities)") + theme_minimal()
print(figS2F)
ggsave(file.path(OUTPUT_DIR, "Figure_S2F_Concat_3Mod.pdf"), figS2F, width = 7, height = 7)
addWorksheet(wb, "FigS2F")
figS2F_data <- data.frame(Embeddings(juxtaData, "concatumap2d"), cellType = juxtaData$layerCellType)
colnames(figS2F_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "FigS2F", figS2F_data)

############### Figure S3A: Spike width vs P2T scatter (colored by cell type) ###############

scatterS3 <- data.frame(
  logWidth = log(as.numeric(juxtaData$width)),
  logP2T = log(as.numeric(juxtaData$ratio_p2t)),
  layerCellType = juxtaData$layerCellType
)
scatterS3 <- scatterS3[is.finite(scatterS3$logWidth) & is.finite(scatterS3$logP2T), ]

figS3A <- ggscatterhist(scatterS3, x = "logWidth", y = "logP2T",
                        color = "layerCellType", margin.params = list(fill = "layerCellType"),
                        main.plot.size = 2, margin.plot.size = 1,
                        ggtheme = theme_minimal(),
                        title = "Spike Width vs P2T (Cell Type)")
print(figS3A)
pdf(file.path(OUTPUT_DIR, "Figure_S3A_Width_P2T_CellType.pdf"), width = 9, height = 7)
print(figS3A)
dev.off()
addWorksheet(wb, "FigS3A")
writeData(wb, "FigS3A", scatterS3[, c("logWidth", "logP2T", "layerCellType")])

############### Figure S3B: Spike width vs P2T scatter (GMM clusters) ###############

gmmInput <- scatterS3[, c("logWidth", "logP2T")]
gmmFit <- Mclust(gmmInput)
scatterS3$GMMcluster <- as.factor(gmmFit$classification)

figS3B <- ggscatterhist(scatterS3, x = "logWidth", y = "logP2T",
                        color = "GMMcluster", margin.params = list(fill = "GMMcluster"),
                        main.plot.size = 2, margin.plot.size = 1,
                        ggtheme = theme_minimal(),
                        title = "Spike Width vs P2T (GMM)")
print(figS3B)
pdf(file.path(OUTPUT_DIR, "Figure_S3B_Width_P2T_GMM.pdf"), width = 9, height = 7)
print(figS3B)
dev.off()
addWorksheet(wb, "FigS3B")
writeData(wb, "FigS3B", scatterS3[, c("logWidth", "logP2T", "GMMcluster")])

############### Figure S3C: Pie charts for GMM clusters ###############

figS3C_collector <- list()
for (gc in sort(unique(scatterS3$GMMcluster))) {
  idx <- which(scatterS3$GMMcluster == gc)
  origCT <- juxtaData$layerCellType[as.numeric(rownames(scatterS3)[idx])]
  ctTable <- table(origCT) / length(origCT)
  ctDF <- data.frame(cellType = names(ctTable), proportion = as.numeric(ctTable))
  ctDF$GMMcluster <- gc
  figS3C_collector[[length(figS3C_collector) + 1]] <- ctDF

  figS3C <- ggplot(ctDF, aes(x = "", y = proportion, fill = cellType)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") + theme_void() +
    ggtitle(paste("GMM Cluster", gc))
  print(figS3C)
  ggsave(file.path(OUTPUT_DIR, paste0("Figure_S3C_GMM_Cluster_", gc, ".pdf")), figS3C, width = 5, height = 5)
}
addWorksheet(wb, "FigS3C")
writeData(wb, "FigS3C", do.call(rbind, figS3C_collector))

############### Figure S4A: Leiden clustering at different resolutions ###############

figS4A_collector <- list()
resolutions <- c(0.5, 1.0, 1.5, 2.0, 2.5)
for (resV in resolutions) {
  juxtaData <- FindClusters(juxtaData, algorithm = ALGORITHM, resolution = resV, verbose = FALSE)
  tmpDF <- data.frame(Embeddings(juxtaData, "wnn.umap"), cluster = juxtaData$seurat_clusters, resolution = resV)
  colnames(tmpDF)[1:2] <- c("UMAP1", "UMAP2")
  figS4A_collector[[length(figS4A_collector) + 1]] <- tmpDF
  p <- DimPlot(juxtaData, reduction = 'wnn.umap', group.by = "seurat_clusters", pt.size = 2) +
    ggtitle(paste("Leiden res =", resV)) + theme_minimal()
  print(p)
}
addWorksheet(wb, "FigS4A")
writeData(wb, "FigS4A", do.call(rbind, figS4A_collector))

############### Figure S4B: MARI vs resolution for different n_neighbors ###############

neighborValues <- c(2, 5, 20, 30, 40)
ariByNeighbor <- list()

for (nn in neighborValues) {
  tryCatch({
    juxtaDataTemp <- FindMultiModalNeighbors(
      juxtaData,
      reduction.list = list("WFpca", "ISIpca"),
      dims.list = list(pcs, pcs),
      modality.weight.name = c("WF.weight.nn", "ISI.weight.nn"),
      k.nn = nn,
      verbose = FALSE
    )
    juxtaDataTemp <- RunUMAP(juxtaDataTemp, nn.name = "weighted.nn",
                             reduction.name = "wnn.umap.nn",
                             reduction.key = "wnnUMAPnn_",
                             seed.use = UMAP.SEED)
    juxtaDataTemp <- FindNeighbors(juxtaDataTemp, reduction = "wnn.umap.nn", dims = 1:2)

    ariVals <- c()
    for (resV in xV) {
      juxtaDataTemp <- FindClusters(juxtaDataTemp, algorithm = 3, resolution = resV, verbose = FALSE)
      ariVals <- c(ariVals, MARI(juxtaDataTemp$seurat_clusters, juxtaDataTemp$CellType))
    }
    ariByNeighbor[[as.character(nn)]] <- ariVals
  }, error = function(e) {
    message(paste("Skipping k.nn =", nn, ":", conditionMessage(e)))
  })
}

if (length(ariByNeighbor) > 0) {
  successNN <- names(ariByNeighbor)
  ariNeighborDF <- data.frame(
    resolution = rep(xV, length(successNN)),
    MARI = unlist(ariByNeighbor),
    n_neighbors = rep(paste0("n=", successNN), each = length(xV))
  )
  figS4B <- ggplot(ariNeighborDF, aes(x = resolution, y = MARI, color = n_neighbors)) +
    geom_line(linewidth = 1) +
    theme_minimal() + ggtitle("MARI vs Resolution (varying n_neighbors)") +
    xlab("Resolution") + ylab("MARI")
  print(figS4B)
} else {
  message("FigS4B: No neighbor values succeeded, skipping.")
}

############### Figure S5A: Marker size on 3-modality PhysMAP ###############

wnnEmbed3 <- Embeddings(juxtaData, reduction = "wnn.umap.3mod")
sizeDF3 <- data.frame(
  UMAP1 = wnnEmbed3[, 1], UMAP2 = wnnEmbed3[, 2],
  width = as.numeric(juxtaData$width),
  ratio_p2t = as.numeric(juxtaData$ratio_p2t),
  latency = as.numeric(juxtaData$latency)
)

figS5A_width <- ggplot(sizeDF3, aes(x = UMAP1, y = UMAP2, size = width)) +
  geom_point(alpha = 0.7) + scale_size_continuous(range = c(0.5, 5)) +
  theme_minimal() + ggtitle("3-Mod PhysMAP: Spike Width")
print(figS5A_width)

figS5A_p2t <- ggplot(sizeDF3, aes(x = UMAP1, y = UMAP2, size = ratio_p2t)) +
  geom_point(alpha = 0.7) + scale_size_continuous(range = c(0.5, 5)) +
  theme_minimal() + ggtitle("3-Mod PhysMAP: P2T Ratio")
print(figS5A_p2t)

sizeDF3$logLatency <- log(sizeDF3$latency + 1)
figS5A_latency <- ggplot(sizeDF3, aes(x = UMAP1, y = UMAP2, size = logLatency)) +
  geom_point(alpha = 0.7) + scale_size_continuous(range = c(0.5, 5)) +
  theme_minimal() + ggtitle("3-Mod PhysMAP: Onset Latency (log)")
print(figS5A_latency)
addWorksheet(wb, "FigS5A")
writeData(wb, "FigS5A", sizeDF3)

############### Figure S5B: 3-modality PhysMAP cell type + Leiden ###############

juxtaData <- FindNeighbors(juxtaData, reduction = "wnn.umap.3mod", dims = 1:2)
juxtaData <- FindClusters(juxtaData, algorithm = ALGORITHM, resolution = RESOLUTION, verbose = FALSE)

figS5B_celltype <- DimPlot(juxtaData, reduction = 'wnn.umap.3mod',
                           group.by = "layerCellType", pt.size = 2) +
  ggtitle("3-Mod PhysMAP: Cell Type") + theme_minimal()
print(figS5B_celltype)

figS5B_clusters <- DimPlot(juxtaData, reduction = 'wnn.umap.3mod',
                           group.by = "seurat_clusters", pt.size = 2) +
  ggtitle("3-Mod PhysMAP: Leiden (res=2)") + theme_minimal()
print(figS5B_clusters)
addWorksheet(wb, "FigS5B_CellType")
figS5B_ct <- data.frame(Embeddings(juxtaData, "wnn.umap.3mod"), cellType = juxtaData$layerCellType)
colnames(figS5B_ct)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "FigS5B_CellType", figS5B_ct)

addWorksheet(wb, "FigS5B_Clusters")
figS5B_cl <- data.frame(Embeddings(juxtaData, "wnn.umap.3mod"), cluster = juxtaData$seurat_clusters)
colnames(figS5B_cl)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "FigS5B_Clusters", figS5B_cl)

############### Figure S5C: MARI scores for 3-modality PhysMAP ###############

pooledARI3 <- c()
for (resV in xV) {
  juxtaData <- FindClusters(juxtaData, algorithm = 3, resolution = resV, verbose = FALSE)
  pooledARI3 <- c(pooledARI3, MARI(juxtaData$seurat_clusters, juxtaData$CellType))
}

ariDF3 <- data.frame(
  resolution = rep(xV, 5),
  MARI = c(isiARI, wfARI, psthARI, concatARI, pooledARI3),
  Modality = rep(c("ISI", "WF", "PSTH", "Concat", "PhysMAP (3-mod)"), each = length(xV))
)

figS5C <- ggplot(ariDF3, aes(x = resolution, y = MARI, color = Modality)) +
  geom_line(linewidth = 1) +
  theme_minimal() + ggtitle("MARI vs Resolution (3 Modalities)") +
  xlab("Resolution") + ylab("MARI")
print(figS5C)
ggsave(file.path(OUTPUT_DIR, "Figure_S5C_MARI_3Mod.pdf"), figS5C, width = 8, height = 6)
addWorksheet(wb, "FigS5C")
writeData(wb, "FigS5C", ariDF3)

############### Figure S5D: GMM MARI vs number of components ###############
# Sweep GMM components for each modality UMAP and PhysMAP WNN UMAP
# Repeated over multiple iterations to compute SEM error bars

compRange <- 2:8
nIterations <- 50

gmmReductions <- list(
  WF       = "WFumap2d",
  ISI      = "ISIumap2d",
  PSTH     = "PSTHumap2d",
  PhysMAP  = "wnn.umap"
)

# Collect all per-iteration MARI values
gmmRawDF <- data.frame()
for (modName in names(gmmReductions)) {
  embedGMM <- Embeddings(juxtaData, reduction = gmmReductions[[modName]])
  embedMat <- embedGMM[, 1:2]
  for (iter in seq_len(nIterations)) {
    set.seed(iter * 101)
    # Use random pairs for EM initialization to get variability across runs
    randPairs <- hcRandomPairs(embedMat)
    gmmARI <- c()
    for (nComp in compRange) {
      gmmRes <- Mclust(embedMat, G = nComp,
                       initialization = list(hcPairs = randPairs))
      gmmARI <- c(gmmARI, MARI(as.factor(gmmRes$classification), juxtaData$CellType))
    }
    gmmRawDF <- rbind(gmmRawDF, data.frame(
      components = compRange,
      MARI = gmmARI,
      Modality = modName,
      iteration = iter
    ))
  }
}

# Compute mean and SEM across iterations
gmmSummaryDF <- gmmRawDF %>%
  dplyr::group_by(components, Modality) %>%
  dplyr::summarise(
    meanMARI = mean(MARI),
    sdMARI   = sd(MARI),
    nIter    = dplyr::n(),
    semMARI  = sd(MARI) / sqrt(dplyr::n()),
    .groups  = "drop"
  )

figS5D <- ggplot(gmmSummaryDF, aes(x = components, y = meanMARI, color = Modality)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = meanMARI - semMARI, ymax = meanMARI + semMARI),
                width = 0.3, linewidth = 0.6) +
  theme_minimal() + ggtitle("GMM Clustering: MARI vs Components") +
  xlab("Number of Components") + ylab("MARI")
print(figS5D)
ggsave(file.path(OUTPUT_DIR, "Figure_S5D_GMM_MARI.pdf"), figS5D, width = 9, height = 6)
addWorksheet(wb, "FigS5D")
writeData(wb, "FigS5D", gmmSummaryDF)

############### Figure S5E: Classification accuracy (all modalities + raw + derived) ###############

# Compute high-dimensional (10D) WNN UMAP for 3-modality classification
juxtaData <- RunUMAP(juxtaData, nn.name = "weighted.nn",
                     reduction.name = "wnn.umap.3mod.hd",
                     reduction.key = "wnnUMAP3HD_",
                     seed.use = UMAP.SEED,
                     n.components = UMAP.components)

wnnEmbed3mod <- data.frame(Embeddings(juxtaData, reduction = "wnn.umap.3mod.hd"))

allResultsS5E <- list()
modNamesS5E <- c("PhysMAP (3-mod)", "WF", "ISI", "PSTH", "Features")
embedListS5E <- list(wnnEmbed3mod, wfEmbed2F, isiEmbed2F, psthEmbed2F, featDF2F)

for (i in seq_along(modNamesS5E)) {
  res <- doClassifyJuxta(embedListS5E[[i]], juxtaData, numreps = 5,
                         method = 'boot', seedV = 1, whichType = 'layercells')
  res$uF$Modality <- modNamesS5E[i]
  allResultsS5E[[i]] <- res$uF
}

classDFS5E <- do.call(rbind, allResultsS5E)
classDFS5E$cellClass <- gsub("Class: ", "", classDFS5E$cellClass)

figS5E <- ggplot(classDFS5E, aes(x = cellClass, y = AccV, color = Modality, group = Modality)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  theme_minimal() + ggtitle("Classification Accuracy (3 Modalities)") +
  xlab("Cell Type") + ylab("Balanced Accuracy (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(figS5E)
ggsave(file.path(OUTPUT_DIR, "Figure_S5E_Classification_3Mod.pdf"), figS5E, width = 9, height = 6)
addWorksheet(wb, "FigS5E")
writeData(wb, "FigS5E", classDFS5E)

############### Figure S6A: Classifier accuracy across embedding dimensions ###############

dimValues <- c(30, 20, 10, 5, 2)
dimResults <- list()

for (nd in dimValues) {
  tryCatch({
    juxtaDataDim <- FindMultiModalNeighbors(
      juxtaData,
      reduction.list = list("WFpca", "ISIpca"),
      dims.list = list(1:nd, 1:nd),
      modality.weight.name = c("WF.weight.dim", "ISI.weight.dim"),
      verbose = FALSE
    )
    juxtaDataDim <- RunUMAP(juxtaDataDim, nn.name = "weighted.nn",
                            reduction.name = "wnn.umap.dim",
                            reduction.key = "wnnUMAPdim_",
                            seed.use = UMAP.SEED,
                            n.components = UMAP.components)
    embedDim <- data.frame(Embeddings(juxtaDataDim, reduction = "wnn.umap.dim"))
    res <- doClassifyJuxta(embedDim, juxtaDataDim, numreps = 5,
                           method = 'boot', seedV = 1, whichType = 'layercells')
    res$uF$nDims <- paste0("d=", nd)
    dimResults[[as.character(nd)]] <- res$uF
  }, error = function(e) {
    message(paste("Skipping nDims =", nd, ":", conditionMessage(e)))
  })
}

if (length(dimResults) > 0) {
  dimDF <- do.call(rbind, dimResults)
  dimDF$cellClass <- gsub("Class: ", "", dimDF$cellClass)
  figS6A <- ggplot(dimDF, aes(x = cellClass, y = AccV, color = nDims, group = nDims)) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    theme_minimal() + ggtitle("Accuracy vs Embedding Dimensions") +
    xlab("Cell Type") + ylab("Balanced Accuracy (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(figS6A)
  ggsave(file.path(OUTPUT_DIR, "Figure_S6A_Accuracy_Dims.pdf"), figS6A, width = 9, height = 6)
  addWorksheet(wb, "FigS6A")
  writeData(wb, "FigS6A", dimDF)
} else {
  message("FigS6A: No dimension values succeeded, skipping.")
}

############### Figure S6B: Classifier accuracy across classifier types ###############

wnnEmbedS6B <- data.frame(Embeddings(juxtaData, reduction = "wnn.umap"))
tempCells <- str_trim(juxtaData$layerCellType)
idx <- tempCells %in% c("E-4", "E-5", "FS-4", "FS-5", "SOM-nan")
trainDF <- wnnEmbedS6B[idx, ]
trainDF$origCells <- factor(tempCells[idx])

set.seed(1)
i <- createDataPartition(trainDF$origCells, times = 1, p = 0.8, list = FALSE)
training <- trainDF[i[, 1], ]
testingset <- trainDF[-i[, 1], ]

classifierMethods <- c("gbm", "rf", "svmRadial", "rpart", "knn", "nnet")
classifierNames <- c("GBM", "RF", "SVM", "CT", "KNN", "NN")
classResults <- list()

for (j in seq_along(classifierMethods)) {
  tryCatch({
    ctrl <- trainControl(method = "boot", number = 5)
    model <- train(origCells ~ ., data = training, method = classifierMethods[j],
                   trControl = ctrl, verbose = FALSE)
    Rpred <- confusionMatrix(predict(model, newdata = testingset), testingset$origCells)
    U <- data.frame(Rpred$byClass)
    U <- U[c(3, 1, 5, 2, 4), ]
    uF <- data.frame(cellClass = rownames(U), AccV = U$Balanced.Accuracy * 100,
                     Classifier = classifierNames[j])
    classResults[[j]] <- uF
  }, error = function(e) {
    message(paste("Skipping classifier", classifierNames[j], ":", conditionMessage(e)))
  })
}

classResults <- Filter(Negate(is.null), classResults)
if (length(classResults) > 0) {
  classTypeDF <- do.call(rbind, classResults)
  classTypeDF$cellClass <- gsub("Class: ", "", classTypeDF$cellClass)
  figS6B <- ggplot(classTypeDF, aes(x = cellClass, y = AccV, color = Classifier, group = Classifier)) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    theme_minimal() + ggtitle("Accuracy vs Classifier Type") +
    xlab("Cell Type") + ylab("Balanced Accuracy (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(figS6B)
  ggsave(file.path(OUTPUT_DIR, "Figure_S6B_Accuracy_Classifiers.pdf"), figS6B, width = 9, height = 6)
  addWorksheet(wb, "FigS6B")
  writeData(wb, "FigS6B", classTypeDF)
} else {
  message("FigS6B: No classifiers succeeded, skipping.")
}

############### Figure S7: PCA variance explained scree plots ###############

wfPCA <- Embeddings(juxtaData, reduction = "WFpca")
isiPCA <- Embeddings(juxtaData, reduction = "ISIpca")
psthPCA <- Embeddings(juxtaData, reduction = "PSTHpca")

wfVar <- Stdev(juxtaData, reduction = "WFpca")^2
isiVar <- Stdev(juxtaData, reduction = "ISIpca")^2
psthVar <- Stdev(juxtaData, reduction = "PSTHpca")^2

screeDF <- data.frame(
  PC = rep(seq_along(wfVar), 3),
  VarExplained = c(wfVar / sum(wfVar) * 100,
                   isiVar / sum(isiVar) * 100,
                   psthVar / sum(psthVar) * 100),
  Modality = rep(c("WF", "ISI", "PSTH"), each = length(wfVar))
)

figS7 <- ggplot(screeDF, aes(x = PC, y = VarExplained, color = Modality)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  theme_minimal() + ggtitle("PCA Variance Explained") +
  xlab("Principal Component") + ylab("% Variance Explained")
print(figS7)
ggsave(file.path(OUTPUT_DIR, "Figure_S7_PCA_Scree.pdf"), figS7, width = 8, height = 5)
addWorksheet(wb, "FigS7")
writeData(wb, "FigS7", screeDF)

############### Save Source Data ###############

saveWorkbook(wb, file.path(SOURCE_DATA_DIR, "01_Figures_1_2_S2-S7_SourceData.xlsx"), overwrite = TRUE)
