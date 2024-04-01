basedir <- dirname(sys.frame(1)$ofile)
setwd(basedir)

source("../constants.R")
source("./helperFunctions.r")

library(ggExtra)
library(ggpubr)
library(scatterpie)
library(reticulate)

pcs = 1:30
nDims = 30
numComponents = UMAP.components #originally 10

# Calculates the individual representations and plots them nicely.
allData = readJianingData("./JianingData/MergedData.mat");
juxtaData =  allData$data;
tempFeat = calcRepresentation(juxtaData, 'features',5,1:5)
juxtaData =  tempFeat$data;
pF = tempFeat$p1

tempISI = calcRepresentation(juxtaData, 'ISI',nDims, pcs, "ISI", TRUE, nc=numComponents, metric=UMAP.metric)
juxtaData =  tempISI$data;
pI = tempISI$p1
ARIvCell = tempISI$layerCellTypeARI;
ARIvCellalone = tempISI$cellTypeARI

tempWF = calcRepresentation(juxtaData, 'WF', nDims, pcs, "WF", FALSE, nc=numComponents, metric=UMAP.metric)
juxtaData =  tempWF$data;
pW = tempWF$p1
ARIvCell = cbind(ARIvCell, tempWF$layerCellTypeARI);
ARIvCellalone = cbind(ARIvCellalone, tempWF$cellTypeARI)

tempPSTH = calcRepresentation(juxtaData, 'PSTH', nDims, pcs, "PSTH",TRUE, nc=numComponents, metric=UMAP.metric)
juxtaData =  tempPSTH$data;
pP = tempPSTH$p1
ARIvCell = cbind(ARIvCell, tempPSTH$layerCellTypeARI);
ARIvCellalone = cbind(ARIvCellalone, tempPSTH$cellTypeARI)

tempConcat = calcRepresentation(juxtaData, 'concat', nDims, pcs, "concat",TRUE, nc=numComponents, metric=UMAP.metric)
juxtaData =  tempConcat$data;
pC = tempConcat$p1
ARIvCell = cbind(ARIvCell, tempConcat$layerCellTypeARI);
ARIvCellalone = cbind(ARIvCellalone, tempConcat$cellTypeARI)

pInd = pW+theme_void() | pI+theme_void() | pP+theme_void() | pC+theme_void()

embed.WF <- as.data.frame(tempWF$data@reductions$WFumap2d@cell.embeddings)
clusts.2res.WF <- as.data.frame(tempWF$data@meta.data$WF_snn_res.2)
colnames(clusts.2res.WF) <- c('clust.id')
clust.embed.WF <- cbind(embed.WF,clusts.2res.WF)
pWF.clust <- ggplot(clust.embed.WF,aes(x=wfumap2d_1,y=wfumap2d_2,group=clust.id,col=clust.id))+geom_point()+theme_void()

embed.ISI <- as.data.frame(tempISI$data@reductions$ISIumap2d@cell.embeddings)
clusts.2res.ISI <- as.data.frame(tempISI$data@meta.data$ISI_snn_res.2)
colnames(clusts.2res.ISI) <- c('clust.id')
clust.embed.ISI <- cbind(embed.ISI,clusts.2res.ISI)
pISI.clust <- ggplot(clust.embed.ISI,aes(x=isiumap2d_1,y=isiumap2d_2,group=clust.id,col=clust.id))+geom_point()+theme_void()

embed.PSTH <- as.data.frame(tempPSTH$data@reductions$PSTHumap2d@cell.embeddings)
clusts.2res.PSTH <- as.data.frame(tempPSTH$data@meta.data$PSTH_snn_res.2)
colnames(clusts.2res.PSTH) <- c('clust.id')
clust.embed.PSTH <- cbind(embed.PSTH,clusts.2res.PSTH)
pPSTH.clust <- ggplot(clust.embed.PSTH,aes(x=psthumap2d_1,y=psthumap2d_2,group=clust.id,col=clust.id))+geom_point()+theme_void()

pComb.clust <- pWF.clust+ggtitle("WF") | pISI.clust+ggtitle("ISI") | pPSTH.clust+ggtitle("PSTH")

# Calculates merged representations.
juxtaData <- FindMultiModalNeighbors(
  juxtaData, reduction.list = list("WFpca","ISIpca","PSTHpca"), 
  dims.list = list(1:30, 1:30,1:30)
)
juxtaData <- RunUMAP(juxtaData, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                   reduction.key = "wnnUMAP_", 
                   seed.use=UMAP.SEED, 
                   metric="euclidean", 
                   n.neighbors = UMAP.neighbors,
                     )

data <- FindClusters(juxtaData,graph.name = "wsnn", algorithm = ALGORITHM, resolution = 2, verbose = FALSE)
p4 <- DimPlot(juxtaData, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 8, pt.size=2)
p4 = p4 + theme_minimal()

xV = seq(0.1,3,0.1);
repeats = seq(1,2);

for (i in repeats){
  d1 = c()
  d2 = c()
  for(resV in xV)
  {
      data <- FindClusters(juxtaData,graph.name = "wsnn", 
                           algorithm = ALGORITHM, 
                           resolution = resV, verbose = FALSE, random.seed=sample(1:100000,1, replace=T))
      d1= c(d1,MARI(data$seurat_clusters, data$CellType))
      d2 = c(d2, MARI(data$seurat_clusters, data$layerCellType))
  }
  if (i == 1){
    d1_runs = d1
    d2_runs = d2
  }
  d1_runs=cbind(d1_runs,d1)
  d2_runs=cbind(d2_runs,d2)
}

ARIvCell = cbind(ARIvCell, d2, xV);
ARIvCellalone = cbind(ARIvCellalone, d1, xV)

library(reshape2)
ARIv = data.frame(ARIvCell)
colnames(ARIv) = c('ISI','WF','PSTH','Concat','Pooled','xV')
ARIvmelted = melt(ARIv, id.vars = "xV")
ARIg1 = ggplot(data=ARIvmelted, aes(x=xV, y=value, group=variable, col=variable)) + geom_line(size=2)
ARIg1 = ARIg1 + theme_minimal()

ARIg1

ARIv = data.frame(ARIvCellalone)
colnames(ARIv) = c('ISI','WF','PSTH','Pooled','xV')
ARIvmelted = melt(ARIv, id.vars = "xV")
ARIg2 = ggplot(data=ARIvmelted, aes(x=xV, y=value, group=variable, col=variable)) + geom_line(size=2)
ARIg2 = ARIg2 + theme_minimal()


p1 <- DimPlot(juxtaData, reduction = 'wnn.umap', 
              group.by = "layerCellType", pt.size=2, repel=TRUE)
p1 = p1 + theme_void()

p2 <- DimPlot(juxtaData, reduction = 'wnn.umap', 
              group.by = "CellType", pt.size=2)
p2 = p2 + theme_void()

juxtaData <- FindClusters(juxtaData,graph.name = "wsnn", algorithm = 2, resolution = RESOLUTION, verbose = FALSE)
p3 <- DimPlot(juxtaData, reduction = 'wnn.umap', label = TRUE, repel = TRUE, pt.size=2)
p3 = p3 + theme_void()


pComb = p1 | p2 | p3

show(pInd)
show(pComb)

load("./JianingData/width.Rda")
load("./JianingData/ratio_p2t.Rda")

# Show plots of points colored in different ways according to other variables
# such as latency
E = Embeddings(juxtaData[["wnn.umap"]])
umapEmbeddings = data.frame(E, allData$F, juxtaData$WF.weight, juxtaData$ISI.weight, juxtaData$PSTH.weight, width=width*1000, ratio_p2t);
  
p = ggplot(umapEmbeddings, aes(x=wnnUMAP_1, y=wnnUMAP_2)) + geom_point(aes(color=layerCellType, size=width)) 
p = p + scale_size_continuous(range=c(0.5,5))
p = p + theme_minimal() + ggtitle("Width of Waveform")
pWidth = p


p = ggplot(umapEmbeddings, aes(x=wnnUMAP_1, y=wnnUMAP_2)) + geom_point(aes(color=layerCellType, size=ratio_p2t))
p = p + scale_size_continuous(range=c(.5,5))
p = p + theme_minimal() + ggtitle("Peak to Trough")
pP2t = p

show(pWidth | pP2t)

pWidth.hist <- ggplot(umapEmbeddings,aes(x=width,fill=layerCellType)) + geom_histogram(alpha=0.5,position="identity") + scale_x_log10() 

umapEmbeddings$log.width = log(umapEmbeddings$width)
umapEmbeddings$log.ratio_p2t = log(umapEmbeddings$ratio_p2t)
ggscatterhist(umapEmbeddings, x="log.width", y= "log.ratio_p2t", 
              color="layerCellType", margin.plot = "histogram",
              margin.params = list(fill = "layerCellType", color = "white", size = 0.3))

umapEmbeddings$region <- factor(1:length(umapEmbeddings$Depth))
p.pie <- ggplot() + geom_scatterpie(aes(x=wnnUMAP_1, y=wnnUMAP_2, group=region), 
                                    cols=c("juxtaData.WF.weight", "juxtaData.ISI.weight", "juxtaData.PSTH.weight"),data=umapEmbeddings)

show(p.pie)

p.hist <- ggplot() + 
  geom_histogram(aes(x=juxtaData.WF.weight), binwidth=0.05, fill="#F39922", alpha=0.6, data = umapEmbeddings) +
  geom_histogram(aes(x=juxtaData.PSTH.weight), binwidth=0.05, fill="#0B78BE", alpha=0.6, data = umapEmbeddings) +
  geom_histogram(aes(x=juxtaData.ISI.weight), binwidth=0.05, fill="#12A84B", alpha=0.6, data = umapEmbeddings) +
  theme_minimal()

show(p.hist)

SOM.p2tratio <- umapEmbeddings$log.ratio_p2t[umapEmbeddings$layerCellType %in% c("SOM-nan")]
SOM.width <- umapEmbeddings$log.width[umapEmbeddings$layerCellType %in% c("SOM-nan")]

umapEmbeddings$log.latency = log(umapEmbeddings$latency)

p = ggplot(umapEmbeddings, aes(x=wnnUMAP_1, y=wnnUMAP_2)) + geom_point(aes(color=layerCellType, size=log.latency))
p = p + scale_size_continuous(range=c(.5,5))
p = p + theme_minimal() + ggtitle("Latency")