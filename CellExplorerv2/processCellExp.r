basedir <- dirname(sys.frame(1)$ofile)
setwd(basedir)

library(tidyverse)
library(caret)
library(nnet)
UMAP.SEED = 3;
ALGORITHM = 3;


rawD = readMat("./Data/theta_modulation_index.mat")
thetaMod = rawD$test

featureData = readMat("../lookupTable/data/features_cellExp_final_Dec2023.mat")
featureData = readMat("./Data/features_cellExp_new_Nov2023.mat")

features = featureData$features

acgData = readMat("./Data/CellExplorer_ACG.mat")

WFce = readMat("./Data/finalWaveforms.mat");
X_waveform = WFce$X
cType = WFce$CellTypeNames

## Now look at area*type labels
area.cType = unlist(featureData$area)
area.cType = sub("VISal", "VIS", area.cType)
area.cType = sub("VISam", "VIS", area.cType)
area.cType = sub("VISl", "VIS", area.cType)
area.cType = sub("VISp", "VIS", area.cType)
area.cType = sub("VISrl", "VIS", area.cType)
area.cType = sub("OTHER", "VIS", area.cType)
area.cType = sub("CA1", "HIP", area.cType)
area.cType = sub("CA3", "HIP", area.cType)
area.cType = data.frame(cbind(matrix(area.cType),cType))
area.cType = area.cType %>% unite("X3", X1:X2, remove = TRUE)

dataSize = dim(X_waveform)
cellIds = c(seq(1, dataSize[1]))
rownames(X_waveform) = cellIds
colnames(X_waveform) = c(seq(1,dataSize[2]))
data = CreateSeuratObject(counts = t(X_waveform), assay = "WF")
data@meta.data <-cbind(data@meta.data,cType)
data@meta.data <-cbind(data@meta.data,area.cType$X3)
#data@meta.data = cbind(data@meta.data, thetaMod)
fS = unlist(featureData$source)
data@meta.data = cbind(data@meta.data, fS)
fSfull = unlist(featureData$sourceFull)
data@meta.data = cbind(data@meta.data, fSfull)

area = unlist(featureData$area)
data@meta.data = cbind(data@meta.data, area)


load('./Data/isi_cellExp.Rda')
X_ISI = t(isi)
dataSize = dim(X_ISI)
rownames(X_ISI) = cellIds
colnames(X_ISI) = c(seq(1,dataSize[2]))
ISI_assay <- CreateAssayObject(counts = t(X_ISI))
data[["ISI1"]] = ISI_assay


acg = acgData$acgW
X_ACG = t(acg)
dataSize = dim(X_ACG)
rownames(X_ACG) = cellIds
colnames(X_ACG) = c(seq(1,dataSize[2]))
ISI_assay <- CreateAssayObject(counts = t(X_ACG))
data[["ACG"]] = ISI_assay


X_features = features
dataSize = dim(X_features)
rownames(X_features) = cellIds
colnames(X_features) = c(seq(1,dataSize[2]))
features_assay <- CreateAssayObject(counts = t(X_features))
data[["features"]] = features_assay

xV = seq(0.05,3,0.05);

numpcs = 40
dimV = 1:40
cValues = c("magenta","red", "orange","red", "blue","cyan","yellow","purple","violet");
RESOLUTION = 1;

# cValues = c("blue", "red","lightgray");

# Finish loading data

runAnalysis = function(data, whichAssay, RESOLUTION=0.75, 
                       metricV = "cosine", normalize=FALSE,
                       numPCS=numpcs, 
                       DimV=dimV, nc=10)
{
  DefaultAssay(data) = whichAssay
  if(normalize)
  {
    data <- NormalizeData(data, normalization.method = "CLR", margin=2)
  }
  data <- ScaleData(data)
  data <- FindVariableFeatures(data)
  data <- RunPCA(data, verbose = FALSE, reduction.name=paste0(whichAssay, 'pca'),npcs=numPCS)
  data <- RunPCA(data, verbose = FALSE, npcs=numPCS)
  data <- FindNeighbors(data, dims=DimV)
  data <- RunUMAP(data, dims=DimV, reduction.name=paste0(whichAssay, 'umap'), 
                  metric=metricV, n.components = nc, seed.use = UMAP.SEED)
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
data = runAnalysis(data, whichAssay, RESOLUTION = RESOLUTION, normalize=TRUE, numPCS=10, DimV=1:10)
p1 = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2)
pF = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2, group.by = "cType", 
                cols=cValues) + theme_minimal() + ggtitle("Features")
show(pF)
dF = calcARI(data, "features")


whichAssay = "ACG"
data = runAnalysis(data, whichAssay, RESOLUTION = RESOLUTION, normalize=TRUE)
p1 = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2)
pACG = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2, group.by = "cType", 
             cols=cValues) + theme_minimal() + ggtitle("ACG")
show(pACG)
dA = calcARI(data, "ACG")



whichAssay = "ISI1"
data = runAnalysis(data, whichAssay, RESOLUTION = RESOLUTION, normalize=TRUE)
p1 = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2)

pISI1 = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2, group.by = "cType", 
              cols=cValues) + theme_minimal() + ggtitle("ISI")
dISI = calcARI(data, "ISI1")


whichAssay = "WF"
data = runAnalysis(data, whichAssay, RESOLUTION = RESOLUTION, normalize=TRUE)
dWF = calcARI(data, whichAssay)

p1 = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2)
pWF = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2, group.by = "cType", 
              cols=cValues) + theme_minimal() + ggtitle("WF")

pComb = pWF | pISI1 | pF | pACG
show(pComb)



data <- FindMultiModalNeighbors(
  data, reduction.list = list("WFpca","ISI1pca","ACGpca","featurespca"), 
  dims.list = list(1:40, 1:40,1:40,1:10)
)
data <- RunUMAP(data, nn.name = "weighted.nn", 
                reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use=3)

data <- FindClusters(data, graph.name = "wsnn", algorithm = 2, resolution = RESOLUTION, verbose = FALSE)

p33 <- DimPlot(data, reduction = 'wnn.umap', pt.size=2, label = FALSE, 
               repel = TRUE, label.size = 8)
p33 <- p33 + theme_minimal()
show(p33)

p32 <- DimPlot(data, reduction = 'wnn.umap', pt.size=2, label = FALSE, 
               repel = TRUE, label.size = 8, group.by="cType", cols = cValues)
p32 <- p32 + theme_minimal()
show(p32)

p34 <- DimPlot(data, reduction = 'wnn.umap', pt.size=2, label = FALSE, 
               repel = TRUE, label.size = 8, group.by="fS", cols = cValues)
p34 <- p34 + theme_minimal()
show(p34)

area.cValues = c("red","red","red","blue","blue","blue","blue","blue","blue")

p35 <- DimPlot(data, reduction = 'wnn.umap', pt.size=2, label = FALSE, 
               repel = TRUE, label.size = 8, group.by="area", cols = area.cValues)
p35 <- p35 + theme_minimal()
show(p35)



dWNN = c()
for(resV in xV)
{
  data <- FindClusters(data,graph.name = "wsnn", algorithm = ALGORITHM, 
                       resolution = resV, verbose = FALSE)
  dWNN= c(dWNN,MARI(data$seurat_clusters, data$cType))
  
}

ARIvCell = cbind(dWF, dISI, dA, dF, dWNN, xV);

library(reshape2)
ARIv = data.frame(ARIvCell)
colnames(ARIv) = c('WF','ISI','ACG','Features','Pooled','xV')
ARIvmelted = melt(ARIv, id.vars = "xV")
ARIg1 = ggplot(data=ARIvmelted, aes(x=xV, y=value, group=variable, col=variable)) + geom_line(size=2)
ARIg1 = ARIg1 + theme_minimal()


## Finish plotting single and also 




plotProb = function(data, whichAssay, modulation, cType)
{
  E = Embeddings(data[[whichAssay]])
  umapEmbeddings = data.frame(E,modulation, cType);
  colnames(umapEmbeddings) = c("UMAP_1","UMAP_2","modulation","cType")
  
  p = ggplot(umapEmbeddings, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(color=cType, size=modulation)) 
  p = p + scale_size_continuous(range=c(0.5,5)) + scale_color_manual(values = cValues) #scale_color_manual(values=c("#000000", "#999999", "#E69F00", "#56B4E9"))
  p = p + theme_minimal()
  pProb = p
  return(pProb)
}

pJoint = plotProb(data, "wnn.umap", t(thetaMod), cType)

p36 <- DimPlot(data, reduction = 'wnn.umap', pt.size=2, label = FALSE, 
               repel = TRUE, label.size = 8, group.by="area.cType$X3")
p36 <- p36 + theme_minimal()
show(p36)