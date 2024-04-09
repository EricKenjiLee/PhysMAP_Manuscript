basedir <- dirname(sys.frame(1)$ofile)
setwd(basedir)

here::i_am("README.md")

readJianingData = function(MatFile)
{

  # Reads a MAT file containing data from Jianing Yu's study on mouse S1
  # Should contain a file called "MergedData.mat"
  
  temp = readMat(MatFile)
  numCells = dim(temp$depth);
  cellIds = c(seq(numCells[2],numCells[1]))
  depthV = temp$depth
  colnames(depthV) = "Depth"
  rownames(depthV) = cellIds
  depth = depthV
  
  latency = temp$latency
  colnames(latency) = "latency"
  rownames(latency) = cellIds
  latency[is.na(latency)] = 0;
  
  onset = temp$onset.rate
  colnames(onset) = "onset";
  rownames(onset) = cellIds;
  onset[is.na(onset)] = 0;
  
  
  mergedISI = temp$ISI;
  rownames(mergedISI) = cellIds
  colnames(mergedISI) = c(seq(1,100))
  data = CreateSeuratObject(counts = t(mergedISI), assay = "ISI")
  data$orig.ident = temp$true.cell.type;
  
  load(here::here("juxtacellular","JianingData","width.Rda"))
  load(here::here("juxtacellular","JianingData","ratio_p2t.Rda"))
  ratio.p2t = ratio_p2t
  temp$features = cbind(temp$features,width,ratio.p2t)

  otherFeatures = temp$features
  otherFeatures[is.na(otherFeatures)] = 0;
  rownames(otherFeatures) = cellIds
  colnames(otherFeatures) = c("Depth","OnsetRate","Latency","StimFR","StimFF","width","ratio.p2t");
  depth_assay = CreateAssayObject(counts = t(otherFeatures), assay = "otherFeatures")
  data[["features"]] = depth_assay
  
  WF = temp$waveformV
  rownames(WF) = cellIds
  colnames(WF) = c(seq(1,351))
  WF_assay <- CreateAssayObject(counts = t(WF))
  data[["WF"]] = WF_assay
  
  PSTH = cbind(temp$AllResp)
  rownames(PSTH) = cellIds
  colnames(PSTH) = c(seq(1,151))
  PSTH_assay <- CreateAssayObject(counts = t(PSTH))
  data[["PSTH"]] = PSTH_assay
  
  concat = cbind(WF, mergedISI, PSTH)
  rownames(concat) = cellIds
  colnames(concat) = c(seq(1,602))
  concat_assay <- CreateAssayObject(counts = t(concat))
  data[["concat"]] = concat_assay
  
  layerData = temp$layer;
  colnames(layerData) = "Layer";
  
  cellType = temp$true.cell.type;
  colnames(cellType) = "CellType"
  layerCellType = paste0(trimws(cellType),'-', trimws(layerData));
  
  data@meta.data <-cbind(data@meta.data,layerData)
  data@meta.data <-cbind(data@meta.data,cellType)
  data@meta.data <-cbind(data@meta.data,depthV)
  data@meta.data <-cbind(data@meta.data,latency)
  data@meta.data <-cbind(data@meta.data,layerCellType)
  data@meta.data <-cbind(data@meta.data, onset);
  data@meta.data <-cbind(data@meta.data, width)
  data@meta.data <-cbind(data@meta.data, ratio_p2t)
  
  F = data.frame(depth = depthV);
  F$latency = latency
  F$onset = onset;
  F$cellType = cellType
  F$layer = layerData
  F$layerCellType = layerCellType
  F$width = width
  F$ratio_p2t = ratio.p2t
  
  width[is.na(width)] = 0
  ratio.p2t[is.na(ratio.p2t)] = 0
  
  otherFeatures = cbind(depth, latency, onset, width, ratio.p2t)
  D = prcomp(x = otherFeatures, scale=TRUE, center=TRUE)
  F$pc1 = D$x[,1]
  F$pc2 = D$x[,2]
  F$pc3 = D$x[,3]
  F$pc4 = D$x[,4]
  F$pc5 = D$x[,5]
  
  
  return(list(data = data, F=F));
  

}

calcRepresentation = function(data, whichAssay, 
                              numpcs=50, 
                              dimV=1:50, 
                              plotTitle="", normalize=TRUE,
                              metric=UMAP.metric, nc=2)
{
  # Does all the heavy lifting -- 
  #
  # 
  
  DefaultAssay(data) = whichAssay
  if(normalize)
      data <- NormalizeData(data,method="CLR", margin = norm.margin)
  data <- FindVariableFeatures(data)
  data <- ScaleData(data)
  data <- RunPCA(data, verbose = FALSE,  
                 reduction.name=paste0(whichAssay, 'pca'),
                 npcs=numpcs)
  data <- RunPCA(data, verbose = FALSE, npcs=numpcs)
  data <- FindNeighbors(data, dims=dimV)
  data <- RunUMAP(data, dims=dimV, reduction.name=paste0(whichAssay, 'umap'),
                  metric = metric, n.components = nc, seed.use = UMAP.SEED)
  data <- RunUMAP(data, dims=dimV, reduction.name=paste0(whichAssay, 'umap2d'), metric = metric)
  data <- RunUMAP(data, dims=dimV,  metric = metric)
  data <- FindClusters(data, algorithm = ALGORITHM, 
                       resolution = RESOLUTION, verbose = FALSE)
  
  p1 = DimPlot(data, reduction = 'umap',  group.by = "layerCellType", pt.size=2) + ggtitle(plotTitle)
  p1 = p1 + theme_minimal()
  
  d1 = c()
  d2 = c()
  xV = seq(0.1,3,0.1);
  for(resV in xV)
  {
    data <- FindClusters(data, algorithm = 3, resolution = resV, verbose = FALSE)
    d1= c(d1,MARI(data$seurat_clusters, data$CellType))
    d2 = c(d2, MARI(data$seurat_clusters, data$layerCellType))
  }
  
  
  ari1 = MARI(data$seurat_clusters, data$CellType)
  ari2 = MARI(data$seurat_clusters, data$layerCellType)
  
  return(list(data = data, p1 = p1, cellTypeARI=d1, layerCellTypeARI = d2));
}


