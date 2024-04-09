library(here)
library(caret)

basedir <- dirname(sys.frame(1)$ofile)
setwd(basedir)

here::i_am("README.md")

source(here::here("constants.R"))

numpcs = 29
dimV = 1:29
cValues = c("orange","blue","red","lightgray");
NC = 5

load(here::here("InvivoA1","A1data","jaramillo_celltypes.Rda"))
load(here::here("InvivoA1","A1data","spikeWaveformsNormalized.Rda"))
M = readMat(here::here("InvivoA1","A1data","santiagoQC.mat"))

load(here::here("InvivoA1","A1data","isiViolations.Rda"))
load(here::here("InvivoA1","A1data","ISI_dist.Rda"))
ISIv = readMat(here::here("InvivoA1","A1data","santiagoISI.mat"))

spikeQuality = M$X1[9,]
selV = spikeQuality >= median(spikeQuality)
selV = spikeQuality >= 4
pValues = (M$X1[4,])/(M$X1[1,]+2)
pValues = -log(M$X1[8,])


X_waveform = spikeWaveformsNormalized[selV ,]
cType = jaramillo_celltypes[selV]
# cType[cType=="EXC"] = "undef"


# X_waveform = spikeWaveformsNormalized;
dataSize = dim(X_waveform)
cellIds = c(seq(1, dataSize[1]))
rownames(X_waveform) = cellIds
colnames(X_waveform) = c(seq(1,dataSize[2]))
data = CreateSeuratObject(counts = t(X_waveform), assay = "WF")

data@meta.data <-cbind(data@meta.data,cType)

X_ISI1 = t(ISI_dist[,selV])
X_ISI2 = ISIv$data[selV ,]

spikeAmplitudes = read.csv(here::here("InvivoA1","A1data","spikeAmplitudes_Fixed.csv"))
spikeAmplitudes = spikeAmplitudes[selV,]

X_ISI = X_ISI1
dataSize = dim(X_ISI)
rownames(X_ISI) = cellIds
colnames(X_ISI) = c(seq(1,dataSize[2]))
ISI_assay <- CreateAssayObject(counts = t(X_ISI))
data[["ISI1"]] = ISI_assay



spikeAmplitudes = read.csv(here::here("InvivoA1","A1data","spikeAmplitudes_Fixed.csv"))
spikeAmplitudes = spikeAmplitudes[selV,]

X_ISI = X_ISI2
dataSize = dim(X_ISI)
rownames(X_ISI) = cellIds
colnames(X_ISI) = c(seq(1,dataSize[2]))
ISI_assay <- CreateAssayObject(counts = t(X_ISI))
data[["ISI2"]] = ISI_assay

# cValues = c("blue", "red","lightgray");

runAnalysis = function(data, whichAssay, RESOLUTION=0.75, 
                       metricV = "cosine", nc=NC)
{
  DefaultAssay(data) = whichAssay
  data <- NormalizeData(data, normalization.method = "CLR", margin=2)
  data <- ScaleData(data)
  data <- FindVariableFeatures(data)
  data <- RunPCA(data, verbose = FALSE, reduction.name=paste0(whichAssay, 'pca'),
                 npcs=numpcs)
  data <- RunPCA(data, verbose = FALSE, npcs=numpcs)
  data <- FindNeighbors(data, dims=dimV)
  data <- RunUMAP(data, dims=dimV, reduction.name=paste0(whichAssay, 'umap'), 
                  metric=metricV, seed.use = UMAP.SEED, 
                  n.components = nc)
  data <- RunUMAP(data, dims=dimV, metric = metricV, 
                  reduction.name=paste0(whichAssay, 'plotumap'))
  data <- FindClusters(data, algorithm = ALGORITHM, resolution = RESOLUTION, verbose = FALSE)
  return(data)
}


xV = seq(0.1,3,0.2);



RESOLUTION = 1
whichAssay = "WF"
data = runAnalysis(data, whichAssay, RESOLUTION = RESOLUTION)
p1 = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2)
pWF = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2, 
              group.by = "cType", cols = cValues)


dWF = c()
iX = data$cType %in% c("PV","SOM","undef")
for(resV in xV)
{
  data <- FindClusters(data, algorithm = ALGORITHM, resolution = resV, verbose = FALSE)
  dWF= c(dWF,MARI(data$seurat_clusters[iX], data$cType[iX]))
}


whichAssay = "ISI1"
data = runAnalysis(data, whichAssay)
p1 = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2)
pISI1 = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2, 
              group.by = "cType", cols = cValues)
dISI = c()
for(resV in xV)
{
  data <- FindClusters(data, algorithm = ALGORITHM, resolution = resV, verbose = FALSE)
  dISI = c(dISI,MARI(data$seurat_clusters[iX], data$cType[iX]))
  
}

whichAssay = "ISI2"
data = runAnalysis(data, whichAssay)
p1 = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2)
pISI1 = DimPlot(data, reduction = paste0(whichAssay, 'plotumap'), pt.size=2, 
                group.by = "cType", cols = cValues)


data <- FindMultiModalNeighbors(
  data, reduction.list = list("WFpca","ISI2pca"), 
  dims.list = list(dimV, dimV, dimV)
)
data <- RunUMAP(data, nn.name = "weighted.nn", 
                reduction.name = "wnn.umap", 
                reduction.key = "wnnUMAP_", seed.use=UMAP.SEED)
data <- FindClusters(data, graph.name = "wsnn", algorithm = ALGORITHM, 
                     resolution = RESOLUTION, verbose = FALSE)
p32 <- DimPlot(data, reduction = 'wnn.umap', label = FALSE, repel = TRUE, 
               label.size = 8, group.by="cType", cols = cValues, pt.size = 2)
p32 + theme_minimal()


p31 <- DimPlot(data, reduction = 'wnn.umap', repel = TRUE, label.size = 8, pt.size = 4)
p31 + theme_minimal()

iX = data$cType %in% c("PV","SOM","undef")
dWNN = c()
for(resV in xV)
{
  data <- FindClusters(data,graph.name = "wsnn", algorithm = ALGORITHM, 
                       resolution = resV, verbose = FALSE)
  dWNN= c(dWNN,MARI(data$seurat_clusters[iX], data$cType[iX]))
  
}

ARIvCell = cbind(dWF, dISI, dWNN, xV);

library(reshape2)
ARIv = data.frame(ARIvCell)
colnames(ARIv) = c('WF','ISI','Pooled','xV')
ARIvmelted = melt(ARIv, id.vars = "xV")
ARIg1 = ggplot(data=ARIvmelted, aes(x=xV, y=value, group=variable, col=variable)) + geom_line(size=2)
ARIg1 = ARIg1 + theme_minimal()


plotProb = function(data, whichAssay, modulation, cType)
{
  E = Embeddings(data[[whichAssay]])
  umapEmbeddings = data.frame(E,modulation, cType);
  colnames(umapEmbeddings) = c("UMAP_1","UMAP_2","modulation","cType")
  
  p = ggplot(umapEmbeddings, aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(color=cType, size=modulation)) 
  p = p + scale_size_continuous(range=c(0.5,5)) + scale_color_manual(values = cValues) #scale_color_manual(values=c("#000000", "#999999", "#E69F00", "#56B4E9"))
  p = p + theme_void()
  pProb = p
  return(pProb)
}

pJoint = plotProb(data, "wnn.umap", pValues[selV], cType)
pWF = plotProb(data, "WFplotumap", pValues[selV], cType)
pISI = plotProb(data, "ISI2plotumap", pValues[selV], cType)
show(pJoint | pWF | pISI)


pJoint = plotProb(data, "wnn.umap", isiViolations[selV], cType)
pWF = plotProb(data, "WFplotumap", isiViolations[selV], cType)
show(pJoint | pWF)


calcAccuracy = function(data, whichUMAP, method, numreps=5, numIter=20, p=0.7)
{
  if(whichUMAP %in% c("wnn.umap","WFumap", "ISI1umap","ISI2umap", "wnn.umap2"))
  {
    print("using UMAP")
    load(here::here("InvivoA1","A1data","spikeWidth.Rda"))
    load(here::here("InvivoA1","A1data","isiViolations.Rda"))
    
    spikeAmplitudes = read.csv(here::here("InvivoA1","A1data","spikeAmplitudes_Fixed.csv"))
    
    spikeAmplitudes = spikeAmplitudes[selV,]
    
    E = data.frame(Embeddings(data[[whichUMAP]]),spikeAmplitudes)
    # E = data.frame(Embeddings(data[[whichUMAP]]))
    print(dim(E))
  }else
  {
    load(here::here("InvivoA1","A1data","spikeWidth.Rda"))
    load(here::here("InvivoA1","A1data","isiViolations.Rda"))
    
    spikeAmplitudes = read.csv(here::here("InvivoA1","A1data","spikeAmplitudes_Fixed.csv"))
    
    spikeAmplitudes = spikeAmplitudes[selV,]
    
    E = data.frame(spikeWidth[selV], spikeAmplitudes)
    print(dim(E))
  }
  
  tempCells = data$cType
  selector = tempCells %in% c('PV','SOM','undef');
  reOrgcells = unlist(factor(tempCells[selector]))
  E = E[selector,]
  E$origCells = unlist(reOrgcells)
  
  for(i in 1:numIter)
  {
    
    method = method
    numreps = numreps;
    p = p;
    set.seed(i)
    
    split <- createDataPartition(E$origCells, times = 1, p = 0.75, list = FALSE)
    training = E[split[,1],]
    testingset = E[-split[,1],]
    

    ctrl <- trainControl(method = method, number=numreps)
    #fit a regression model and use k-fold CV to evaluate performance
    model <- train(origCells~., data = E, method = "gbm", 
                   trControl = ctrl, verbose=FALSE)
    print(mean(model$results$Accuracy))
    
    Rpred = confusionMatrix(predict(model, newdata = testingset), testingset$origCells)
    E1 = data.frame(Rpred$byClass)
    Rpred$byClass
    # browser()
  
    accV = E1$Balanced.Accuracy;
    print(accV)
    
    if(i==1)
      AccVall = accV
    else
      AccVall = rbind(AccVall, accV)
    
  }
  return(AccVall)
}

data <- RunUMAP(data, nn.name = "weighted.nn", 
                   reduction.name = "wnn.umap2", 
                   reduction.key = "wnnUMAP2_", seed.use=UMAP.SEED, n.components = NC,
)



AccComb = calcAccuracy(data, "wnn.umap2", "repeatedcv", numIter=50)
AccWF = calcAccuracy(data, "WFumap", "repeatedcv", numIter=50)
AccFeatures = calcAccuracy(data, "features", "repeatedcv", numIter=50)
AccISI = calcAccuracy(data, "ISI1umap", "repeatedcv", numIter=50)


# AccWF = calcAccuracy(data, "WFumap", "cv", numIter=10)


dataSummary <- function(data, varname, groupnames){
  require(plyr)
  summaryFunc <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  dataSum<-ddply(data, groupnames, .fun=summaryFunc,
                 varname)
  dataSum <- rename(dataSum, c("mean" = varname))
  return(dataSum)
}


nreps = dim(AccWF)[1]
rawAccData = AccWF
colnames(rawAccData) = c("PV","SOM","UNDEF")
rownames(rawAccData) = seq(1,nreps)
rawAccDataWF = melt(rawAccData)
rawAccDataWF["Type"] = "WF"

rawAccData = AccComb
colnames(rawAccData) = c("PV","SOM","UNDEF")
rownames(rawAccData) = seq(1,nreps)
rawAccDataComb = melt(rawAccData)
rawAccDataComb["Type"] = "WNN"

rawAccData = AccISI
colnames(rawAccData) = c("PV","SOM","UNDEF")
rownames(rawAccData) = seq(1,nreps)
rawAccDataISI = melt(rawAccData)
rawAccDataISI["Type"] = "ISI"

rawAccData = AccFeatures
colnames(rawAccData) = c("PV","SOM","UNDEF")
rownames(rawAccData) = seq(1,nreps)
rawAccDataFeatures = melt(rawAccData)
rawAccDataFeatures["Type"] = "Features"

combData = rbind(rawAccDataWF, rawAccDataComb, rawAccDataISI, rawAccDataFeatures)

colnames(combData) = c("Run","CellType","Acc","Modality")

library(ggthemes)
summaryData = dataSummary(combData, varname="Acc", groupnames = c("CellType","Modality"))
summaryData$se = summaryData$sd/sqrt(nreps)
p<- ggplot(summaryData, aes(x=CellType, y=Acc, group=Modality, color=Modality)) + 
  geom_line() + theme_classic() +  
  geom_point(aes(size=8)) + theme(text = element_text(size=20)) + ylim(0.7,1.0) + 
  geom_errorbar(aes(ymin=Acc-se, ymax=Acc+se), width=.2,
                position=position_dodge(0.0))

show(p)
