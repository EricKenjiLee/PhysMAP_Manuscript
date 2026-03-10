# classifyDataFigure.r
# Two-panel classifier performance figure for the PhysMAP manuscript.
#
# Top panel:  GBM classifier across 7 representations at their native
#             embedding dimensionality (PhysMAP at UMAP.components,
#             individual UMAP modalities at 30, Features at 2D, etc.)
#
# Bottom panel: 6 classifiers (GBM, RF, SVM, CART, KNN, Neural Net)
#               across the same 7 representations, all at embedding dim 30.
#
# Prerequisite: processJianing.r must have been sourced first so that
#               `juxtaData` exists in the environment with WF/ISI/PSTH
#               UMAP reductions already computed.

library(tidyverse)
library(caret)
library(nnet)
library(gbm)
library(randomForest)
library(kernlab)
library(rpart)
library(reshape2)
library(cowplot)

set.seed(42)

basedir <- dirname(sys.frame(1)$ofile)
setwd(basedir)

here::i_am("README.md")

source(here::here("constants.R"))

# ================================================================
# Helper Functions
# ================================================================

dataSummary <- function(data, varname, groupnames){
  summaryFunc <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  dataSum <- plyr::ddply(data, groupnames, .fun=summaryFunc, varname)
  dataSum <- plyr::rename(dataSum, c("mean" = varname))
  return(dataSum)
}

doClassify = function(E, seuratDat, numreps=5,
                      method='repeatedcv', seedV=1, repeats=10,
                      whichType='layercells', classifierMethod='gbm')
{
  if(whichType == 'layercells') {
    tempCells = str_trim(seuratDat$layerCellType)
    idx = tempCells %in% c("E-4","E-5","FS-4","FS-5","SOM-nan")
    E = E[idx,]
    tempCells = tempCells[idx]
    origCells = factor(tempCells)
  } else {
    origCells = factor(str_trim(seuratDat$CellType))
    tempCells = str_trim(seuratDat$CellType)
    idx = tempCells %in% c("E","FS","SOM")
    E = E[idx,]
    tempCells = tempCells[idx]
    origCells = factor(tempCells)
  }
  
  E$origCells = origCells
  set.seed(seedV)
  i <- createDataPartition(E$origCells, times=1, p=0.8, list=FALSE)
  training  = E[i[,1],]
  testingset = E[-i[,1],]
  
  ctrl <- trainControl(method=method, number=numreps)
  
  if(classifierMethod == 'gbm') {
    model <- train(origCells~., data=training, method=classifierMethod,
                   trControl=ctrl, verbose=FALSE)
  } else if(classifierMethod == 'nnet') {
    model <- train(origCells~., data=training, method=classifierMethod,
                   trControl=ctrl, trace=FALSE, MaxNWts=5000)
  } else {
    model <- train(origCells~., data=training, method=classifierMethod,
                   trControl=ctrl)
  }
  
  Rpred = confusionMatrix(predict(model, newdata=testingset), testingset$origCells)
  Acc = Rpred$overall[1]
  
  U = data.frame(Rpred$byClass)
  U = U[c(3,1,5,2,4),]
  uF = data.frame(cellClass = rownames(U), AccV = U$Balanced.Accuracy*100)
  
  return(list(uF=uF, AccV=Acc))
}

collectResults = function(accMatrix, cellClasses, nreps, modalityName) {
  rawAccData = t(accMatrix)
  colnames(rawAccData) = cellClasses
  rownames(rawAccData) = seq(1, nreps)
  df = melt(rawAccData)
  df["Type"] = modalityName
  return(df)
}

cleanCellType = function(x) {
  x = gsub("^Class: ", "", x)
  x = gsub("SOM-nan", "SOM", x)
  return(x)
}

# ================================================================
# Prepare Embeddings
# ================================================================

load('./JianingData/width.Rda')
load('./JianingData/ratio_p2t.Rda')

# PCA of concatenated raw data (WF + ISI + PSTH)
concatRaw = cbind(t(GetAssayData(juxtaData, assay="WF")),
                  t(GetAssayData(juxtaData, assay="ISI")),
                  t(GetAssayData(juxtaData, assay="PSTH")))
variances = apply(concatRaw, 2, var)
concatRaw = concatRaw[, variances > 0]
concatPCA = prcomp(concatRaw, center=TRUE, scale.=TRUE)

# WNN UMAP at native dim (UMAP.components) for top panel
juxtaData <- RunUMAP(juxtaData, nn.name="weighted.nn",
                     reduction.name="wnn.umap.native",
                     reduction.key="wnnUMAPnat_",
                     seed.use=UMAP.SEED,
                     n.components=UMAP.components,
                     metric="correlation")

# WNN UMAP at dim 30 for bottom panel
juxtaData <- RunUMAP(juxtaData, nn.name="weighted.nn",
                     reduction.name="wnn.umap30",
                     reduction.key="wnnUMAP30_",
                     seed.use=UMAP.SEED,
                     n.components=30,
                     metric="correlation")

# ================================================================
# Configuration
# ================================================================

representationNames = c("PhysMAP", "WF", "ISI", "PSTH",
                        "Waveform Metric", "Concatenated", "Raw Data")

classifierMethods = c("gbm", "rf", "svmRadial", "rpart", "knn", "nnet")
classifierLabels  = c("gbm"="GBM", "rf"="Random Forest",
                      "svmRadial"="Radial Basis SVM",
                      "rpart"="Classification Tree",
                      "knn"="KNN", "nnet"="Neural Network")

nIterations = 20

getEmbeddings = function(representation, dimMode='native') {
  if (representation == "PhysMAP") {
    if (dimMode == 'native') {
      E = data.frame(Embeddings(juxtaData[["wnn.umap.native"]]))
    } else {
      E = data.frame(Embeddings(juxtaData[["wnn.umap30"]]))
    }
    E[is.na(E)] = 0
  } else if (representation == "WF") {
    E = data.frame(Embeddings(juxtaData[["WFumap"]]))
  } else if (representation == "ISI") {
    E = data.frame(Embeddings(juxtaData[["ISIumap"]]))
  } else if (representation == "PSTH") {
    E = data.frame(Embeddings(juxtaData[["PSTHumap"]]))
  } else if (representation == "Waveform Metric") {
    E = data.frame(ratio_p2t, width)
    E[is.na(E)] = 0
  } else if (representation == "Concatenated") {
    if (dimMode == 'native') {
      E = data.frame(concatPCA$x[, 1:UMAP.components])
    } else {
      E = data.frame(concatPCA$x[, 1:30])
    }
  } else if (representation == "Raw Data") {
    if (dimMode == 'native') {
      E = data.frame(concatPCA$x[, 1:nrow(concatRaw)])
    } else {
      E = data.frame(concatPCA$x[, 1:30])
    }
  }
  return(E)
}

# ================================================================
# TOP PANEL: GBM at native embedding dimensions
# ================================================================

print("=== Top Panel: GBM at native embedding dimensions ===")

allResults_top = list()

for(rep in representationNames) {
  print(paste("  Representation:", rep))
  accMatrix = NULL
  for(i in 1:nIterations) {
    E = getEmbeddings(rep, dimMode='native')
    currAcc = doClassify(E, juxtaData, seedV=i, classifierMethod='gbm')
    if(is.null(accMatrix)) {
      accMatrix = currAcc$uF$AccV
    } else {
      accMatrix = cbind(accMatrix, currAcc$uF$AccV)
    }
  }
  allResults_top[[rep]] = collectResults(accMatrix, currAcc$uF$cellClass,
                                         nIterations, rep)
}

combData_top = do.call(rbind, allResults_top)
colnames(combData_top) = c("Run", "CellType", "Acc", "Modality")
combData_top$CellType = cleanCellType(combData_top$CellType)

summaryData_top = dataSummary(combData_top, varname="Acc",
                              groupnames=c("CellType", "Modality"))
summaryData_top$se = summaryData_top$sd / sqrt(nIterations)

p_top <- ggplot(summaryData_top,
                aes(x=CellType, y=Acc, group=Modality, color=Modality)) +
  theme_classic() +
  coord_cartesian(clip="off") +
  geom_point(size=3, position=position_dodge(0.5)) +
  geom_line(linewidth=0.5, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=Acc-se, ymax=Acc+se), width=0.2,
                position=position_dodge(0.5)) +
  theme(text=element_text(size=14)) +
  ylim(50, 100) +
  labs(x="Cell Type", y="Balanced Accuracy (%)",
       title=paste0("GBM (Embedding-D = ", as.character(UMAP.components), ")"))

# ================================================================
# BOTTOM PANEL: Multiple classifiers at embedding dim 30
# ================================================================

print("=== Bottom Panel: Multiple classifiers at embedding dim 30 ===")

allResults_bottom = list()

for(clf in classifierMethods) {
  print(paste("  Classifier:", classifierLabels[clf]))
  for(rep in representationNames) {
    print(paste("    Representation:", rep))
    accMatrix = NULL
    for(i in 1:nIterations) {
      E = getEmbeddings(rep, dimMode='dim30')
      currAcc = doClassify(E, juxtaData, seedV=i, classifierMethod=clf)
      if(is.null(accMatrix)) {
        accMatrix = currAcc$uF$AccV
      } else {
        accMatrix = cbind(accMatrix, currAcc$uF$AccV)
      }
    }
    key = paste0(clf, "_", rep)
    df = collectResults(accMatrix, currAcc$uF$cellClass, nIterations, rep)
    df$Classifier = classifierLabels[clf]
    allResults_bottom[[key]] = df
  }
}

combData_bottom = do.call(rbind, allResults_bottom)
colnames(combData_bottom) = c("Run", "CellType", "Acc", "Modality", "Classifier")
combData_bottom$CellType = cleanCellType(combData_bottom$CellType)

summaryData_bottom = dataSummary(combData_bottom, varname="Acc",
                                 groupnames=c("CellType", "Modality", "Classifier"))
summaryData_bottom$se = summaryData_bottom$sd / sqrt(nIterations)

p_bottom <- ggplot(summaryData_bottom,
                   aes(x=CellType, y=Acc, group=Modality, color=Modality)) +
  theme_classic() +
  coord_cartesian(clip="off") +
  geom_point(size=2, position=position_dodge(0.5)) +
  geom_line(linewidth=0.4, position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin=Acc-se, ymax=Acc+se), width=0.2,
                position=position_dodge(0.5)) +
  theme(text=element_text(size=11),
        strip.text=element_text(size=12, face="bold")) +
  ylim(50, 100) +
  labs(x="Cell Type", y="Balanced Accuracy (%)",
       title="Classifiers at Embedding Dimension = 30") +
  facet_wrap(~ Classifier, ncol=3)

# ================================================================
# Combine and Save
# ================================================================

combined_plot = plot_grid(p_top, p_bottom, ncol=1, rel_heights=c(1, 2),
                          labels=c("A", "B"), label_size=16)

ggsave("classifierFigure.pdf", combined_plot, width=14, height=16, dpi=300)
ggsave("classifierFigure.png", combined_plot, width=14, height=16, dpi=300)

print("Figure saved as classifierFigure.pdf and classifierFigure.png")