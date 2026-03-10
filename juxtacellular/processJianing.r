library(here)
library(ggExtra)
library(ggpubr)
library(scatterpie)
library(reticulate)
library(mclust)

set.seed(42)


library(here)
here::i_am("README.md")

source(here::here("constants.R"))
source(here::here("juxtacellular","helperFunctions.R"))

pcs = 1:30
nDims = 30
numComponents = 30 #originally 10
k.nn = 20

# Calculates the individual representations and plots them nicely.
allData = readJianingData(here::here("juxtacellular","JianingData","MergedData.mat"));
juxtaData =  allData$data;
juxtaData$isExcitatory <- ifelse(grepl("^E", trimws(juxtaData$CellType)), "Excitatory", "Inhibitory")

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
clust.embed.WF$isExcitatory <- tempWF$data$isExcitatory
pWF.clust <- ggplot(clust.embed.WF,aes(x=wfumap2d_1,y=wfumap2d_2,group=clust.id,col=clust.id,shape=isExcitatory))+geom_point()+scale_shape_manual(values=c("Excitatory"=16,"Inhibitory"=16))+theme_void()

embed.ISI <- as.data.frame(tempISI$data@reductions$ISIumap2d@cell.embeddings)
clusts.2res.ISI <- as.data.frame(tempISI$data@meta.data$ISI_snn_res.2)
colnames(clusts.2res.ISI) <- c('clust.id')
clust.embed.ISI <- cbind(embed.ISI,clusts.2res.ISI)
clust.embed.ISI$isExcitatory <- tempISI$data$isExcitatory
pISI.clust <- ggplot(clust.embed.ISI,aes(x=isiumap2d_1,y=isiumap2d_2,group=clust.id,col=clust.id,shape=isExcitatory))+geom_point()+scale_shape_manual(values=c("Excitatory"=16,"Inhibitory"=16))+theme_void()

embed.PSTH <- as.data.frame(tempPSTH$data@reductions$PSTHumap2d@cell.embeddings)
clusts.2res.PSTH <- as.data.frame(tempPSTH$data@meta.data$PSTH_snn_res.2)
colnames(clusts.2res.PSTH) <- c('clust.id')
clust.embed.PSTH <- cbind(embed.PSTH,clusts.2res.PSTH)
clust.embed.PSTH$isExcitatory <- tempPSTH$data$isExcitatory
pPSTH.clust <- ggplot(clust.embed.PSTH,aes(x=psthumap2d_1,y=psthumap2d_2,group=clust.id,col=clust.id,shape=isExcitatory))+geom_point()+scale_shape_manual(values=c("Excitatory"=16,"Inhibitory"=16))+theme_void()

pComb.clust <- pWF.clust+ggtitle("WF") | pISI.clust+ggtitle("ISI") | pPSTH.clust+ggtitle("PSTH")

# Calculates merged representations.
juxtaData <- FindMultiModalNeighbors(
  juxtaData, reduction.list = list("WFpca","ISIpca","PSTHpca"), 
  dims.list = list(1:30, 1:30, 1:30)
)
juxtaData <- RunUMAP(juxtaData, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                   reduction.key = "wnnUMAP_", 
                   seed.use=UMAP.SEED, 
                   metric="euclidean", 
                   n.neighbors = UMAP.neighbors,
                     )

data <- FindClusters(juxtaData,graph.name = "wsnn", algorithm = ALGORITHM, resolution = 2, verbose = FALSE)
wnnE4 = as.data.frame(Embeddings(juxtaData, reduction = 'wnn.umap'))
wnnE4$cluster = juxtaData$seurat_clusters
wnnE4$isExcitatory = juxtaData$isExcitatory
p4 <- ggplot(wnnE4, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = cluster, shape = isExcitatory)) +
  geom_point(size = 2) +
  scale_shape_manual(values = c("Excitatory" = 16, "Inhibitory" = 16))
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


wnnE = as.data.frame(Embeddings(juxtaData, reduction = 'wnn.umap'))
wnnE$layerCellType = juxtaData$layerCellType
wnnE$CellType = juxtaData$CellType
wnnE$isExcitatory = juxtaData$isExcitatory

p1 <- ggplot(wnnE, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = layerCellType, shape = isExcitatory)) +
  geom_point(size = 2) +
  scale_shape_manual(values = c("Excitatory" = 16, "Inhibitory" = 16))
p1 = p1 + theme_void()

p2 <- ggplot(wnnE, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = CellType, shape = isExcitatory)) +
  geom_point(size = 2) +
  scale_shape_manual(values = c("Excitatory" = 16, "Inhibitory" = 16))
p2 = p2 + theme_void()

juxtaData <- FindClusters(juxtaData,graph.name = "wsnn", algorithm = 2, resolution = RESOLUTION, verbose = FALSE)
wnnE$cluster = juxtaData$seurat_clusters
p3 <- ggplot(wnnE, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = cluster, shape = isExcitatory)) +
  geom_point(size = 2) +
  scale_shape_manual(values = c("Excitatory" = 16, "Inhibitory" = 16))
p3 = p3 + theme_void()


pComb = p1 | p2 | p3

show(pInd)
show(pComb)

load(here::here("juxtacellular","JianingData","width.Rda"))
load(here::here("juxtacellular","JianingData","ratio_p2t.Rda"))

# Show plots of points colored in different ways according to other variables
# such as latency
E = Embeddings(juxtaData[["wnn.umap"]])
umapEmbeddings = data.frame(E, allData$F, juxtaData$WF.weight, juxtaData$ISI.weight, juxtaData$PSTH.weight, width=width*1000, ratio_p2t);
umapEmbeddings$isExcitatory = juxtaData$isExcitatory

p = ggplot(umapEmbeddings, aes(x=wnnUMAP_1, y=wnnUMAP_2)) + geom_point(aes(color=layerCellType, size=width, shape=isExcitatory))
p = p + scale_size_continuous(range=c(0.5,5))
p = p + scale_shape_manual(values=c("Excitatory"=16,"Inhibitory"=16))
p = p + theme_minimal() + ggtitle("Width of Waveform")
pWidth = p


p = ggplot(umapEmbeddings, aes(x=wnnUMAP_1, y=wnnUMAP_2)) + geom_point(aes(color=layerCellType, size=ratio_p2t, shape=isExcitatory))
p = p + scale_size_continuous(range=c(.5,5))
p = p + scale_shape_manual(values=c("Excitatory"=16,"Inhibitory"=16))
p = p + theme_minimal() + ggtitle("Peak to Trough")
pP2t = p

p = ggplot(umapEmbeddings, aes(x=wnnUMAP_1, y=wnnUMAP_2)) + geom_point(aes(color=layerCellType, size=latency, shape=isExcitatory))
p = p + scale_size_continuous(range=c(.5,5))
p = p + scale_shape_manual(values=c("Excitatory"=16,"Inhibitory"=16))
p = p + theme_minimal() + ggtitle("Onset Latency (s)")
pLatency = p

show(pWidth | pP2t | pLatency)

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

p = ggplot(umapEmbeddings, aes(x=wnnUMAP_1, y=wnnUMAP_2)) + geom_point(aes(color=layerCellType, size=log.latency, shape=isExcitatory))
p = p + scale_size_continuous(range=c(.5,5))
p = p + scale_shape_manual(values=c("Excitatory"=16,"Inhibitory"=16))
p = p + theme_minimal() + ggtitle("Latency")

df <- data.frame(id = 1:246,celltype = juxtaData@meta.data$layerCellType, cluster_ix = juxtaData@meta.data$wsnn_res.2)
class_counts <- df %>%
  dplyr::group_by(cluster_ix, celltype) %>%
  dplyr::summarise(Count = dplyr::n(), .groups = "drop")

class_counts <- class_counts %>%
  dplyr::group_by(cluster_ix) %>%
  dplyr::mutate(Proportion = Count / sum(Count))

ggplot(class_counts, aes(x = "", y = Proportion, fill = celltype)) +
  geom_bar(stat = "identity", width = 1) +  # Bar chart for pie slices
  coord_polar(theta = "y") +               # Convert to pie chart
  facet_wrap(~cluster_ix) +           # One pie chart per secondary class
  labs(title = "Composition of each cluster by cell type",
       fill = "Primary Class") +
  theme_void() +                           # Clean up the chart
  theme(legend.position = "bottom")

# --- Modality weight pie charts per cluster ---
weightDF = data.frame(
  cluster = juxtaData@meta.data$wsnn_res.2,
  WF = juxtaData$WF.weight,
  ISI = juxtaData$ISI.weight,
  PSTH = juxtaData$PSTH.weight
)

meanWeights = weightDF %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(WF = mean(WF), ISI = mean(ISI), PSTH = mean(PSTH), .groups = "drop")

# Relabel clusters as 1-indexed
meanWeights$cluster = factor(as.integer(as.character(meanWeights$cluster)) + 1)

meanWeightsLong = reshape2::melt(meanWeights, id.vars = "cluster",
                                  variable.name = "Modality",
                                  value.name = "Weight")

pWeightPie = ggplot(meanWeightsLong, aes(x = "", y = Weight, fill = Modality)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~cluster) +
  scale_fill_manual(values = c("WF" = "#F39922", "ISI" = "#12A84B", "PSTH" = "#0B78BE")) +
  labs(title = "Mean modality weight per cluster",
       fill = "Modality") +
  theme_void() +
  theme(legend.position = "bottom",
        text = element_text(size = 14))

show(pWeightPie)

############### Figure S5D: GMM MARI vs number of components ###############
# Sweep GMM components for each modality UMAP and PhysMAP WNN UMAP
# Repeated over multiple iterations to compute SEM error bars
compRange <- 2:10
nIterations <- 1 #Use 50 to replicate publication (it just takes too long when needing to run the full notebook)
gmmReductions <- list(
  WF       = "WFumap2d",
  ISI      = "ISIumap2d",
  PSTH     = "PSTHumap2d",
  Concat   = "concatumap2d",
  PhysMAP  = "wnn.umap"
)
# Collect all per-iteration MARI values
gmmRawDF <- data.frame()
for (modName in names(gmmReductions)) {
  embedGMM <- Embeddings(juxtaData, reduction = gmmReductions[[modName]])
  embedMat <- embedGMM[, 1:2]
  for (iter in seq_len(nIterations)) {
    set.seed(iter * UMAP.SEED)
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
  geom_vline(xintercept = c(2, 4, 6, 8, 10), linetype = "dashed", 
             color = "grey50", alpha = 0.5) +
  annotate("text", x = c(2, 4, 6, 8, 10), y = Inf, 
           label = c("2", "4", "6", "8", "10"),
           vjust = -0.5, size = 3, color = "grey40") +
  theme_minimal() + ggtitle("GMM Clustering: MARI vs Components") +
  xlab("Number of Components") + ylab("MARI")
print(figS5D)

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
figS3A$sp <- figS3A$sp + theme(legend.position = "none")
figS3A$xplot <- figS3A$xplot + theme(legend.position = "none")
figS3A$yplot <- figS3A$yplot + theme(legend.position = "none")
print(figS3A)

############### Figure S3B: Spike width vs P2T scatter (GMM clusters) ###############
gmmInput <- scatterS3[, c("logWidth", "logP2T")]
gmmFit <- Mclust(gmmInput, G=3)
scatterS3$GMMcluster <- as.factor(gmmFit$classification)
figS3B <- ggscatterhist(scatterS3, x = "logWidth", y = "logP2T",
                        color = "GMMcluster", margin.params = list(fill = "GMMcluster"),
                        main.plot.size = 2, margin.plot.size = 1,
                        ggtheme = theme_minimal(),
                        title = "Spike Width vs P2T (GMM)")
figS3B$sp <- figS3B$sp + theme(legend.position = "none")
figS3B$xplot <- figS3B$xplot + theme(legend.position = "none")
figS3B$yplot <- figS3B$yplot + theme(legend.position = "none")
print(figS3B)

############### Figure S3C: Pie charts for GMM clusters ###############
for (gc in sort(unique(scatterS3$GMMcluster))) {
  idx <- which(scatterS3$GMMcluster == gc)
  origCT <- juxtaData$layerCellType[as.numeric(rownames(scatterS3)[idx])]
  ctTable <- table(origCT) / length(origCT)
  ctDF <- data.frame(cellType = names(ctTable), proportion = as.numeric(ctTable))
  figS3C <- ggplot(ctDF, aes(x = "", y = proportion, fill = cellType)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") + theme_void() +
    ggtitle(paste("GMM Cluster", gc))
  print(figS3C)
  write_xlsx(ctDF, paste0("pMetricsClust_gc", gc, ".xlsx"))
}