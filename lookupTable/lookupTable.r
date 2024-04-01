
maxPCs = 3;
dims = 1:maxPCs


cValues = c("orange","blue", "gray","red","magenta","cyan","yellow");

refdata <- FindMultiModalNeighbors(
  refdata, reduction.list = list("WFpca","ISI1pca", 'featurespca'), 
  dims.list = list(1:40, 1:40, 1:8), return.intermediate = TRUE,  knn.range = 50
)
refdata = RunSPCA(refdata, npcs=maxPCs, graph = "wsnn")

refdata.anchors = FindTransferAnchors(reference = refdata, query = mapdata,  
                                      reference.reduction = "spca", dims=dims)

refdata <- RunUMAP(refdata, nn.name = "weighted.nn", 
                reduction.name = "wnn.umap2", reduction.key = "wnnUMAP2_", seed.use=3, return.model = TRUE)
refdata.query <- MapQuery(anchorset = refdata.anchors, reference = refdata, 
                           query = mapdata,
                           refdata = list(celltype = "cType"), 
                           reference.reduction = "spca", reduction.model = "wnn.umap2")
p1 <- DimPlot(refdata, reduction = "wnn.umap2", group.by = "cType", label = TRUE, label.size = 3,
              repel = TRUE, cols=cValues) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(refdata.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE, cols=cValues) + NoLegend() + ggtitle("Query transferred labels")

p1 | p2

iX = (mapdata$cType %in% c('PV','SST','Pyra','Juxt'))
iX = (mapdata$cType %in% c('PV','Pyra','SST'))

table((refdata.query$predicted.celltype[iX] == mapdata$cType[iX]) )

# iX = (mapdata$cType %in% c('Pyra'))
# iX = (mapdata$cType %in% c('PV','SST','VIP'))
# iX = (mapdata$cType %in% c('SST'))

# table(mapdata$cType == 'PV')