library(tidyverse)
library(caret)
library(nnet)
library(reshape2)
library(here)

basedir <- dirname(sys.frame(1)$ofile)
setwd(basedir)
here::i_am("README.md")

source(here("constants.R"))

NORMALIZE = TRUE
set.seed(UMAP.SEED)

rawD = readMat(here("lookupTable","data","theta_modulation_index.mat"))
thetaMod = rawD$test

featureData = readMat(here("lookupTable","data","features_cellExp_final_Dec2023.mat"))
features = featureData$features

nData = dim(features)

# selectedData = featureData$nonAllenNeurons
selectedData = readMat(here("CellExplorerv2","Data","validNeurons.mat"))$allValid.ix

actualData = seq(1, nData[1])
fullAllenData = intersect(actualData, selectedData)

aD = length(fullAllenData)
selectedData = fullAllenData[sample(1:aD[1], size=350, replace=F)]
remData = setdiff(fullAllenData, selectedData)

# selectedData = fullAllenData
# remData = selectedData

WFce = readMat(here("lookupTable","data","finalWaveforms.mat"));

WFce$CellTypeNames = gsub("Juxt","Pyra",WFce$CellTypeNames)
WFce$CellTypeNames = gsub("Axo_","PV",WFce$CellTypeNames)

X_waveform = WFce$X[selectedData,]
cType = WFce$CellTypeNames[selectedData]
thetaMod = thetaMod[selectedData]

# X_waveform = spikeWaveformsNormalized;
dataSize = dim(X_waveform)
cellIds = c(seq(1, dataSize[1]))
rownames(X_waveform) = cellIds
colnames(X_waveform) = c(seq(1,dataSize[2]))
refdata = CreateSeuratObject(counts = t(X_waveform), assay = "WF")
refdata@meta.data <-cbind(refdata@meta.data,cType)
refdata@meta.data = cbind(refdata@meta.data, thetaMod)

load(here("lookupTable","data","isi_cellExp.Rda"))
X_ISI = t(isi)
X_ISI = X_ISI[selectedData, ]
dataSize = dim(X_ISI)
rownames(X_ISI) = cellIds
colnames(X_ISI) = c(seq(1,dataSize[2]))
ISI_assay <- CreateAssayObject(counts = t(X_ISI))
refdata[["ISI1"]] = ISI_assay


X_features = features
dataSize = dim(X_features)
X_features = X_features[selectedData, ]
rownames(X_features) = cellIds
colnames(X_features) = c(seq(1,dataSize[2]))
features_assay <- CreateAssayObject(counts = t(X_features))
refdata[["features"]] = features_assay


#
WFce = readMat(here("lookupTable","data","finalWaveforms.mat"));
WFce$CellTypeNames = gsub("Juxt","Pyra",WFce$CellTypeNames)
WFce$CellTypeNames = gsub("Axo_","PV",WFce$CellTypeNames)

X_waveform = WFce$X[remData,]
cType = WFce$CellTypeNames[remData]
thetaMod = thetaMod[remData]

# X_waveform = spikeWaveformsNormalized;
dataSize = dim(X_waveform)
cellIds = c(seq(1, dataSize[1]))
rownames(X_waveform) = cellIds
colnames(X_waveform) = c(seq(1,dataSize[2]))
mapdata = CreateSeuratObject(counts = t(X_waveform), assay = "WF")
mapdata@meta.data <-cbind(mapdata@meta.data,cType)
mapdata@meta.data = cbind(mapdata@meta.data, thetaMod)

load(here("lookupTable","data","isi_cellExp.Rda"))
X_ISI = t(isi)
X_ISI = X_ISI[remData, ]
dataSize = dim(X_ISI)
rownames(X_ISI) = cellIds
colnames(X_ISI) = c(seq(1,dataSize[2]))
ISI_assay <- CreateAssayObject(counts = t(X_ISI))
mapdata[["ISI1"]] = ISI_assay


X_features = features
dataSize = dim(X_features)
X_features = X_features[remData, ]
rownames(X_features) = cellIds
colnames(X_features) = c(seq(1,dataSize[2]))
features_assay <- CreateAssayObject(counts = t(X_features))
mapdata[["features"]] = features_assay


xV = seq(0.1,3,0.1);

numpcs = 40
dimV = 1:40
cValues = c("orange","blue", "gray","red","magenta","cyan","yellow");
RESOLUTION = 1;

# cValues = c("blue", "red","lightgray");

# Finish loading data

runAnalysis = function(data, whichAssay, RESOLUTION=0.75, 
                       metricV = "cosine", normalize=TRUE,
                       numPCS=numpcs, 
                       DimV=dimV, n.comp = 8)
{
  DefaultAssay(data) = whichAssay
  if(normalize)
  {
    data <- NormalizeData(data, normalization.method = "CLR", margin=1)
  }
  data <- ScaleData(data)
  data <- FindVariableFeatures(data)
  data <- RunPCA(data, verbose = FALSE, reduction.name=paste0(whichAssay, 'pca'),npcs=numPCS)
  data <- RunPCA(data, verbose = FALSE, npcs=numPCS)
  data <- FindNeighbors(data, dims=DimV)
  data <- RunUMAP(data, dims=DimV, reduction.name=paste0(whichAssay, 'umap'), 
                  metric=metricV, n.components = n.comp, seed.use = UMAP.SEED)
  data <- RunUMAP(data, dims=DimV, reduction.name=paste0(whichAssay, 'plotumap'), metric=metricV)
  
  data <- FindClusters(data, algorithm = ALGORITHM, resolution = RESOLUTION, verbose = FALSE)
  return(data)
}

calcARI = function(data, whichAssay, res = xV)
{
  d = c()
  DefaultAssay(data) = whichAssay
  for(resV in res)
  {
    data <- FindClusters(data, algorithm = ALGORITHM, 
                         resolution = resV, verbose = FALSE)
    d= c(d,MARI(data$seurat_clusters, data$cType))
  }
  return(d)
}

whichAssay = "features"
refdata = runAnalysis(refdata, whichAssay, RESOLUTION = RESOLUTION, normalize=NORMALIZE, numPCS=10, DimV=1:10)
p1 = DimPlot(refdata, reduction = paste0(whichAssay, 'plotumap'), pt.size=2)
pF = DimPlot(refdata, reduction = paste0(whichAssay, 'plotumap'), pt.size=2, group.by = "cType", 
                cols=cValues) + theme_minimal() + ggtitle("Features")
show(pF)
dF = calcARI(refdata, "features")


whichAssay = "ISI1"
refdata = runAnalysis(refdata, whichAssay, RESOLUTION = RESOLUTION, normalize=NORMALIZE)
p1 = DimPlot(refdata, reduction = paste0(whichAssay, 'plotumap'), pt.size=2)

pISI1 = DimPlot(refdata, reduction = paste0(whichAssay, 'plotumap'), pt.size=2, group.by = "cType", 
              cols=cValues) + theme_minimal() + ggtitle("ISI")
dISI = calcARI(refdata, "ISI1")


whichAssay = "WF"
refdata = runAnalysis(refdata, whichAssay, RESOLUTION = RESOLUTION, normalize=NORMALIZE)
dWF = calcARI(refdata, whichAssay)

p1 = DimPlot(refdata, reduction = paste0(whichAssay, 'plotumap'), pt.size=2)
pWF = DimPlot(refdata, reduction = paste0(whichAssay, 'plotumap'), pt.size=2, group.by = "cType", 
              cols=cValues) + theme_minimal() + ggtitle("WF")

pComb = pWF | pISI1 | pF
show(pComb)



refdata <- FindMultiModalNeighbors(
  refdata, reduction.list = list("WFpca","ISI1pca","featurespca"), 
  dims.list = list(1:40, 1:40,1:8), knn.range = 100
)
refdata <- RunUMAP(refdata, nn.name = "weighted.nn", 
                reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use=UMAP.SEED)

refdata <- FindClusters(refdata, graph.name = "wsnn", algorithm = ALGORITHM, resolution = RESOLUTION, verbose = FALSE)

p33 <- DimPlot(refdata, reduction = 'wnn.umap', pt.size=2, label = FALSE, repel = TRUE, label.size = 8)
p33 <- p33 + theme_minimal()
show(p33)

p32 <- DimPlot(refdata, reduction = 'wnn.umap', pt.size=2, label = FALSE, repel = TRUE, label.size = 8, group.by="cType", cols = cValues)
p32 <- p32 + theme_minimal()
show(p32)


dWNN = c()
for(resV in xV)
{
  refdata <- FindClusters(refdata,graph.name = "wsnn", algorithm = ALGORITHM, 
                       resolution = resV, verbose = FALSE)
  dWNN= c(dWNN,MARI(refdata$seurat_clusters, refdata$cType))
  
}

ARIvCell = cbind(dWF, dISI, dF, dWNN, xV);

library(reshape2)
ARIv = data.frame(ARIvCell)
colnames(ARIv) = c('WF','ISI','Features','Pooled','xV')
ARIvmelted = melt(ARIv, id.vars = "xV")
ARIg1 = ggplot(data=ARIvmelted, aes(x=xV, y=value, group=variable, col=variable)) + geom_line(size=2)
ARIg1 = ARIg1 + theme_minimal()


## Plot classification probabilities


# plotProb = function(data, whichAssay, modulation, cType)
# {
#   E = Embeddings(data[[whichAssay]])
#   umapEmbeddings = data.frame(E,modulation, cType);
#   colnames(umapEmbeddings) = c("UMAP_1","UMAP_2","modulation","cType")
#   
#   p = ggplot(umapEmbeddings, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(color=cType, size=modulation)) 
#   p = p + scale_size_continuous(range=c(0.5,5)) + scale_color_manual(values = cValues) #scale_color_manual(values=c("#000000", "#999999", "#E69F00", "#56B4E9"))
#   p = p + theme_minimal()
#   pProb = p
#   return(pProb)
# }
# 
# pJoint = plotProb(refdata, "wnn.umap", t(thetaMod), cType)
