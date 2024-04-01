library(caret)
set.seed(42)

temp = readMat("~/Documents/GitHub/PhysMAP_Chand/juxtacellular/JianingData/MergedData.mat")
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

load("/Users/kenjilee/Documents/GitHub/PhysMAP_Chand/juxtacellular/JianingData/width.Rda")
load("/Users/kenjilee/Documents/GitHub/PhysMAP_Chand/juxtacellular/JianingData/ratio_p2t.Rda")
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
E = data.frame(concat)

layerData = temp$layer;
colnames(layerData) = "Layer";

cellType = temp$true.cell.type;
colnames(cellType) = "CellType"

layerCellType = paste0(trimws(cellType),'-', trimws(layerData));

tempCells = str_trim(layerCellType)
idx = tempCells %in% c("E-4","E-5","FS-4", "FS-5", "SOM-nan")

E = E[idx,]
tempCells = tempCells[idx]
origCells = factor(tempCells)
E$origCells = origCells

allAcc = list()
meanAccRaw = list()

for(k in 1:20){
  print(k)
  i <- createDataPartition(E$origCells, times=1, p = 0.8, list = FALSE)
  training = E[i[,1],]
  testingset = E[-i[,1],]
  
  ctrl <- trainControl(method = "boot", number=5)
  #fit a regression model and use k-fold CV to evaluate performance
  model <- train(origCells~., data = training, method = "gbm", 
                 trControl = ctrl, verbose=FALSE)
  mean(model$results$Accuracy)
  predict(model, newdata = testingset)
  
  Rpred = confusionMatrix(predict(model, newdata = testingset), testingset$origCells)
  Acc= Rpred$overall[1]
  
  U = data.frame(Rpred$byClass)
  U = U[c(3,1,5,2,4),]
  uF = data.frame(cellClass = rownames(U), AccV = U$Balanced.Accuracy*100)
  
  currAcc <- list(uF = uF, AccV = Acc)
  
  if (k==1){
    allAcc <-  currAcc$uF$AccV
  }
  
  else{
    allAcc <- cbind(allAcc, currAcc$uF$AccV)
  }

  meanAccRaw[k] = unlist(currAcc$AccV)
}

row_mean <- rowMeans(allAcc)
row_ste <- apply(allAcc, 1, sd, na.rm=TRUE)/sqrt(nreps)
rawData <- data.frame(row_mean)
rawData$cell_type <- c("FS-4","E-4","SOM","E-5","FS-5")
rawData$row_ste <- row_ste

p<- ggplot(rawData, aes(x=cell_type, y=row_mean)) + theme_classic() + 
  geom_point(aes(size=6)) + geom_line(aes(size=0.2)) + 
  geom_errorbar(aes(ymin=row_mean-row_ste, ymax=row_mean+row_ste), width=.2) +
  ylim(50,100)
p