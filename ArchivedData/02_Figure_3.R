############### Figure 3: PhysMAP Applied to A1 Extracellular Data ###############
#
# Dataset: Auditory cortex A1 extracellular recordings (Lakunina et al.)
# n = 373 neurons (PV+, SOM+, Excitatory/untagged)
# Produces: Figures 3A-F
#
##############################################################################
DATA_ROOT <- "~/GitHub/PhysMAP_Manuscript"
OUTPUT_DIR <- "~/GitHub/PhysMAP_Manuscript/InvivoA1"

library(caret)
library(reshape2)
library(ggthemes)
library(openxlsx)

source("~/GitHub/PhysMAP_Manuscript/constants.R")
source("~/GitHub/PhysMAP_Manuscript/juxtacellular/helperFunctions.R")

SOURCE_DATA_DIR <- file.path("~/Desktop/PhysMAP_Figures/output/source_data/")
dir.create(SOURCE_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
wb <- createWorkbook()

DATA_DIR <- file.path(DATA_ROOT, "InvivoA1", "A1data")

numpcs <- 29
dimV <- 1:29
NC <- 5

COLOR_PV <- "#2CA69F"
COLOR_SOM <- "#8B5FC7"
COLOR_EXC <- "#E8870E"
COLOR_UNDEF <- "lightgray"
cValues <- c(COLOR_EXC, COLOR_PV, COLOR_SOM, COLOR_UNDEF)

############### Load Data ###############

load(file.path(DATA_DIR, "jaramillo_celltypes.Rda"))
load(file.path(DATA_DIR, "spikeWaveformsNormalized.Rda"))
M <- readMat(file.path(DATA_DIR, "santiagoQC.mat"))

load(file.path(DATA_DIR, "isiViolations.Rda"))
load(file.path(DATA_DIR, "ISI_dist.Rda"))
ISIv <- readMat(file.path(DATA_DIR, "santiagoISI.mat"))

spikeQuality <- M$X1[9, ]
selV <- spikeQuality >= 4
pValues <- (M$X1[4, ]) / (M$X1[1, ] + 2)

############### Build Seurat Object ###############

X_waveform <- spikeWaveformsNormalized[selV, ]
cType <- jaramillo_celltypes[selV]

dataSize <- dim(X_waveform)
cellIds <- seq(1, dataSize[1])
rownames(X_waveform) <- cellIds
colnames(X_waveform) <- seq(1, dataSize[2])
data <- CreateSeuratObject(counts = t(X_waveform), assay = "WF")
data@meta.data <- cbind(data@meta.data, cType)

X_ISI1 <- t(ISI_dist[, selV])
X_ISI2 <- ISIv$data[selV, ]

rownames(X_ISI1) <- cellIds
colnames(X_ISI1) <- seq(1, dim(X_ISI1)[2])
ISI1_assay <- CreateAssayObject(counts = t(X_ISI1))
data[["ISI1"]] <- ISI1_assay

rownames(X_ISI2) <- cellIds
colnames(X_ISI2) <- seq(1, dim(X_ISI2)[2])
ISI2_assay <- CreateAssayObject(counts = t(X_ISI2))
data[["ISI2"]] <- ISI2_assay

spikeAmplitudes <- read.csv(file.path(DATA_DIR, "spikeAmplitudes_Fixed.csv"))
spikeAmplitudes <- spikeAmplitudes[selV, ]

############### Per-Modality Seurat Processing ###############

runAnalysis <- function(data, whichAssay, RESOLUTION = 0.75,
                        metricV = "cosine", nc = NC) {
  DefaultAssay(data) <- whichAssay
  data <- NormalizeData(data, normalization.method = "CLR", margin = 2)
  data <- ScaleData(data)
  data <- FindVariableFeatures(data)
  data <- RunPCA(data, verbose = FALSE,
                 reduction.name = paste0(whichAssay, "pca"), npcs = numpcs)
  data <- RunPCA(data, verbose = FALSE, npcs = numpcs)
  data <- FindNeighbors(data, dims = dimV)
  data <- RunUMAP(data, dims = dimV,
                  reduction.name = paste0(whichAssay, "umap"),
                  metric = metricV, seed.use = UMAP.SEED, n.components = nc)
  data <- RunUMAP(data, dims = dimV, metric = metricV,
                  reduction.name = paste0(whichAssay, "plotumap"))
  data <- FindClusters(data, algorithm = ALGORITHM,
                       resolution = RESOLUTION, verbose = FALSE)
  return(data)
}

RESOLUTION <- 1

data <- runAnalysis(data, "WF", RESOLUTION = RESOLUTION)
data <- runAnalysis(data, "ISI1")
data <- runAnalysis(data, "ISI2")

############### Weighted Nearest Neighbor Integration ###############

data <- FindMultiModalNeighbors(
  data, reduction.list = list("WFpca", "ISI2pca"),
  dims.list = list(dimV, dimV, dimV)
)

data <- RunUMAP(data, nn.name = "weighted.nn",
                reduction.name = "wnn.umap",
                reduction.key = "wnnUMAP_", seed.use = UMAP.SEED)

data <- FindClusters(data, graph.name = "wsnn", algorithm = ALGORITHM,
                     resolution = RESOLUTION, verbose = FALSE)

data <- RunUMAP(data, nn.name = "weighted.nn",
                reduction.name = "wnn.umap2",
                reduction.key = "wnnUMAP2_", seed.use = UMAP.SEED,
                n.components = NC)

############### Plotting Helper ###############

plotProb <- function(data, whichAssay, modulation, cType) {
  E <- Embeddings(data[[whichAssay]])
  umapEmbeddings <- data.frame(E, modulation, cType)
  colnames(umapEmbeddings) <- c("UMAP_1", "UMAP_2", "modulation", "cType")

  p <- ggplot(umapEmbeddings, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = cType, size = modulation)) +
    scale_size_continuous(range = c(0.5, 5)) +
    scale_color_manual(values = cValues) +
    theme_void()
  return(p)
}

############### Figure 3A: PhysMAP A1 Extracellular ###############

pJoint <- plotProb(data, "wnn.umap", pValues[selV], cType)
show(pJoint)
ggsave(file.path(OUTPUT_DIR, "Figure_3A_PhysMAP_A1.pdf"), pJoint, width = 7, height = 7)

addWorksheet(wb, "Fig3A")
fig3A_data <- data.frame(Embeddings(data[["wnn.umap"]]), pValue = pValues[selV], cellType = cType)
colnames(fig3A_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "Fig3A", fig3A_data)

############### Figure 3B: Classification Accuracy by Modality ###############

calcAccuracy <- function(data, whichUMAP, method, numreps = 5,
                         numIter = 20, p = 0.7) {
  if (whichUMAP %in% c("wnn.umap", "WFumap", "ISI1umap", "ISI2umap", "wnn.umap2")) {
    load(file.path(DATA_DIR, "spikeWidth.Rda"))
    load(file.path(DATA_DIR, "isiViolations.Rda"))
    spikeAmps <- read.csv(file.path(DATA_DIR, "spikeAmplitudes_Fixed.csv"))
    spikeAmps <- spikeAmps[selV, ]
    E <- data.frame(Embeddings(data[[whichUMAP]]), spikeAmps)
  } else {
    load(file.path(DATA_DIR, "spikeWidth.Rda"))
    load(file.path(DATA_DIR, "isiViolations.Rda"))
    spikeAmps <- read.csv(file.path(DATA_DIR, "spikeAmplitudes_Fixed.csv"))
    spikeAmps <- spikeAmps[selV, ]
    E <- data.frame(spikeWidth[selV], spikeAmps)
  }

  tempCells <- data$cType
  selector <- tempCells %in% c("PV", "SOM", "undef")
  reOrgcells <- unlist(factor(tempCells[selector]))
  E <- E[selector, ]
  E$origCells <- unlist(reOrgcells)

  for (i in 1:numIter) {
    set.seed(i)
    split <- createDataPartition(E$origCells, times = 1, p = 0.75, list = FALSE)
    training <- E[split[, 1], ]
    testingset <- E[-split[, 1], ]

    ctrl <- trainControl(method = method, number = numreps)
    model <- train(origCells ~ ., data = E, method = "gbm",
                   trControl = ctrl, verbose = FALSE)

    Rpred <- confusionMatrix(predict(model, newdata = testingset),
                             testingset$origCells)
    E1 <- data.frame(Rpred$byClass)
    accV <- E1$Balanced.Accuracy

    if (i == 1) {
      AccVall <- accV
    } else {
      AccVall <- rbind(AccVall, accV)
    }
  }
  return(AccVall)
}

AccComb <- calcAccuracy(data, "wnn.umap2", "repeatedcv", numIter = 50)
AccWF <- calcAccuracy(data, "WFumap", "repeatedcv", numIter = 50)
AccFeatures <- calcAccuracy(data, "features", "repeatedcv", numIter = 50)
AccISI <- calcAccuracy(data, "ISI1umap", "repeatedcv", numIter = 50)

nreps <- dim(AccWF)[1]

rawAccData <- AccWF
colnames(rawAccData) <- c("PV", "SOM", "UNDEF")
rownames(rawAccData) <- seq(1, nreps)
rawAccDataWF <- melt(rawAccData)
rawAccDataWF["Type"] <- "WF"

rawAccData <- AccComb
colnames(rawAccData) <- c("PV", "SOM", "UNDEF")
rownames(rawAccData) <- seq(1, nreps)
rawAccDataComb <- melt(rawAccData)
rawAccDataComb["Type"] <- "WNN"

rawAccData <- AccISI
colnames(rawAccData) <- c("PV", "SOM", "UNDEF")
rownames(rawAccData) <- seq(1, nreps)
rawAccDataISI <- melt(rawAccData)
rawAccDataISI["Type"] <- "ISI"

rawAccData <- AccFeatures
colnames(rawAccData) <- c("PV", "SOM", "UNDEF")
rownames(rawAccData) <- seq(1, nreps)
rawAccDataFeatures <- melt(rawAccData)
rawAccDataFeatures["Type"] <- "Features"

combData <- rbind(rawAccDataWF, rawAccDataComb, rawAccDataISI, rawAccDataFeatures)
colnames(combData) <- c("Run", "CellType", "Acc", "Modality")

summaryData <- dataSummary(combData, varname = "Acc",
                           groupnames = c("CellType", "Modality"))
summaryData$se <- summaryData$sd / sqrt(nreps)

pClassification <- ggplot(summaryData,
                          aes(x = CellType, y = Acc,
                              group = Modality, color = Modality)) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Acc - se, ymax = Acc + se), width = 0.2) +
  ylim(0.7, 1.0) +
  theme_classic() +
  theme(text = element_text(size = 20))

show(pClassification)
ggsave(file.path(OUTPUT_DIR, "Figure_3B_Classification.pdf"), pClassification, width = 9, height = 6)

addWorksheet(wb, "Fig3B")
writeData(wb, "Fig3B", summaryData)

############### Figure 3C: Waveform UMAP ###############

pWF <- plotProb(data, "WFplotumap", pValues[selV], cType)
show(pWF)
ggsave(file.path(OUTPUT_DIR, "Figure_3C_WF_UMAP.pdf"), pWF, width = 7, height = 7)

addWorksheet(wb, "Fig3C")
fig3C_data <- data.frame(Embeddings(data[["WFplotumap"]]), pValue = pValues[selV], cellType = cType)
colnames(fig3C_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "Fig3C", fig3C_data)

############### Figure 3D: ISI Distribution UMAP ###############

pISI <- plotProb(data, "ISI2plotumap", pValues[selV], cType)
show(pISI)
ggsave(file.path(OUTPUT_DIR, "Figure_3D_ISI_UMAP.pdf"), pISI, width = 7, height = 7)

addWorksheet(wb, "Fig3D")
fig3D_data <- data.frame(Embeddings(data[["ISI2plotumap"]]), pValue = pValues[selV], cellType = cType)
colnames(fig3D_data)[1:2] <- c("UMAP1", "UMAP2")
writeData(wb, "Fig3D", fig3D_data)

############### Figure 3E: Spike Width KDE by Cell Type ###############

load(file.path(DATA_DIR, "spikeWidth.Rda"))
spikeWidthSel <- spikeWidth[selV]
cTypeSel <- cType

widthDF <- data.frame(spikeWidth = spikeWidthSel, cType = cTypeSel)
widthTagged <- widthDF[widthDF$cType %in% c("PV", "SOM", "EXC"), ]

pKDE <- ggplot(widthTagged, aes(x = spikeWidth, fill = cType, color = cType)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 alpha = 0.3, position = "identity") +
  geom_density(linewidth = 1.2, alpha = 0) +
  scale_fill_manual(values = c("EXC" = COLOR_EXC,
                                "PV" = COLOR_PV,
                                "SOM" = COLOR_SOM)) +
  scale_color_manual(values = c("EXC" = COLOR_EXC,
                                 "PV" = COLOR_PV,
                                 "SOM" = COLOR_SOM)) +
  labs(x = "Spike Width (ms)", y = "Density") +
  theme_classic() +
  theme(text = element_text(size = 16))

show(pKDE)
ggsave(file.path(OUTPUT_DIR, "Figure_3E_SpikeWidth_KDE.pdf"), pKDE, width = 7, height = 5)

addWorksheet(wb, "Fig3E")
writeData(wb, "Fig3E", widthTagged)

############### Figure 3F: Spike Width Histogram (All Cells) ###############

pHistAll <- ggplot(widthDF, aes(x = spikeWidth, fill = cType)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("EXC" = COLOR_EXC,
                                "PV" = COLOR_PV,
                                "SOM" = COLOR_SOM,
                                "undef" = COLOR_UNDEF)) +
  labs(x = "Spike Width (ms)", y = "Count") +
  theme_classic() +
  theme(text = element_text(size = 16))

show(pHistAll)
ggsave(file.path(OUTPUT_DIR, "Figure_3F_SpikeWidth_Histogram.pdf"), pHistAll, width = 7, height = 5)

addWorksheet(wb, "Fig3F")
writeData(wb, "Fig3F", widthDF)

############### Composite Figure 3 ###############

show(pJoint | pClassification)
show(pWF | pISI)
show(pKDE | pHistAll)

############### Save Source Data ###############

saveWorkbook(wb, file.path(SOURCE_DATA_DIR, "02_Figure_3_SourceData.xlsx"), overwrite = TRUE)
