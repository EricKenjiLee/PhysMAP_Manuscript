library(plyr)
library(dplyr)

library(tidyverse)
library(caret)
library(nnet)


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

doClassify = function(E, seuratDat, numreps=5, 
                      method='boot', seedV=1, whichType='layercells')
{
 
  if(whichType == 'layercells')
  {
    tempCells = str_trim(seuratDat$layerCellType)
    idx = tempCells %in% c("E-4","E-5","FS-4", "FS-5", "SOM-nan")
    # idx = tempCells %in% c("E4","FS4","SOMnan")
    
    E = E[idx,]
    tempCells = tempCells[idx]
    origCells = factor(tempCells);
  }
  else{
    
  origCells = factor(str_trim(seuratDat$CellType));
  
  tempCells = str_trim(seuratDat$CellType)
  idx = tempCells %in% c("E", "FS", "SOM")
  E = E[idx,]
  tempCells = tempCells[idx]
  origCells = factor(tempCells);
  }
  
  E$origCells = origCells

  # E$origCells = origCells
  set.seed(seedV)
  i <- createDataPartition(E$origCells, times=1, p = 0.8, list = FALSE)
  training = E[i[,1],]
  testingset = E[-i[,1],]

  ctrl <- trainControl(method = method, number=numreps)
  #fit a regression model and use k-fold CV to evaluate performance
  model <- train(origCells~., data = training, method = "gbm", 
                 trControl = ctrl, verbose=FALSE)
  mean(model$results$Accuracy)
  predict(model, newdata = testingset)
  
  Rpred = confusionMatrix(predict(model, newdata = testingset), testingset$origCells)
  Acc= Rpred$overall[1]

  U = data.frame(Rpred$byClass)
  U = U[c(3,1,5,2,4),]
  Rpred
  
  uF = data.frame(cellClass = rownames(U), AccV = U$Balanced.Accuracy*100)
  
  return(list(uF = uF, AccV = Acc))
}


juxtaData <- RunUMAP(juxtaData, nn.name = "weighted.nn", 
                reduction.name = "wnn.umap2", 
                reduction.key = "wnnUMAP2_", seed.use=UMAP.SEED, 
                n.components = UMAP.components, metric = "correlation")


pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                     max = 20, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = ".") 


meanAccWnn = list()
meanAccWf = list()
meanAccConcat = list()

print("Running classification analysis ....")
load('/Users/kenjilee/Documents/GitHub/PhysMAP_Chand/juxtacellular/JianingData/width.Rda')
load('/Users/kenjilee/Documents/GitHub/PhysMAP_Chand/juxtacellular/JianingData/ratio_p2t.Rda')

for(i in 1:20)
{
  # print(i)
 
  
  E = data.frame(Embeddings(juxtaData[["wnn.umap2"]]))
  E[is.na(E)] = 0;
  
  currAccWnn = doClassify(E, juxtaData, seedV = i);
  
  
  E = data.frame(width, ratio_p2t)
  # E = data.frame(Embeddings(juxtaData[["WFumap"]]))
  E[is.na(E)] = 0;
  
  E = data.frame(Embeddings(juxtaData[["WFumap"]]))
  currAccWf = doClassify(E, juxtaData, seedV = i);
  
  E = data.frame(Embeddings(juxtaData[["ISIumap"]]));
  currAccISI = doClassify(E, juxtaData, seedV = i);
  
  E = data.frame(Embeddings(juxtaData[["PSTHumap"]]));
  currAccPSTH = doClassify(E, juxtaData, seedV = i);
  
  E = data.frame(ratio_p2t, width)
  E[is.na(E)] = 0;
  currAccFeatures = doClassify(E, juxtaData, seedV = i);
  
  
  
                 
  if(i==1)
  {
    allWnn = currAccWnn$uF$AccV
    allWf = currAccWf$uF$AccV
    allISI = currAccISI$uF$AccV
    allPSTH = currAccPSTH$uF$AccV
    allFeatures = currAccFeatures$uF$AccV
  }
  else
  {
    # print(currAccWnn$AccV)
    allWnn = cbind(allWnn, currAccWnn$uF$AccV)
    allWf = cbind(allWf, currAccWf$uF$AccV)
    allISI =  cbind(allISI, currAccISI$uF$AccV)
    allPSTH = cbind(allPSTH, currAccPSTH$uF$AccV)
    allFeatures = cbind(allFeatures, currAccFeatures$uF$AccV)
  }
  meanAccWnn[i] = unlist(currAccWnn$AccV)
  meanAccWf[i] = unlist(currAccWf$AccV)
  setTxtProgressBar(pb,i)

}

averageClass = data.frame(currAccWnn$uF$cellClass)
averageClass$wnn = rowMeans(allWnn)
averageClass$wf = rowMeans(allWf)
colnames(averageClass) = c('cellClass','WNN','WF')
acc1 = ggplot(melt(averageClass, id.vars = "cellClass")) + 
  geom_line(aes(x=cellClass, y=value, group=variable, col=variable)) + 
  geom_point(aes(x=cellClass, y=value, group=variable, col=variable)) + theme_minimal()

show(acc1)


# wilcox.test(allWf[1,], allWnn[1,], 
#        alternative = "two.sided", paired = TRUE)

nreps = dim(allWnn)[2]
rawAccData = t(allWnn)
colnames(rawAccData) = currAccWnn$uF$cellClass
rownames(rawAccData) = seq(1,nreps)
rawAccDataWNN = melt(rawAccData)
rawAccDataWNN["Type"] = "WNN"

rawAccData = t(allWf)
colnames(rawAccData) = currAccWnn$uF$cellClass
rownames(rawAccData) = seq(1,nreps)
rawAccDataWF = melt(rawAccData)
rawAccDataWF["Type"] = "WF"

rawAccData = t(allISI)
colnames(rawAccData) = currAccWnn$uF$cellClass
rownames(rawAccData) = seq(1,nreps)
rawAccDataISI = melt(rawAccData)
rawAccDataISI["Type"] = "ISI"

rawAccData = t(allPSTH)
colnames(rawAccData) = currAccWnn$uF$cellClass
rownames(rawAccData) = seq(1,nreps)
rawAccDataPSTH = melt(rawAccData)
rawAccDataPSTH["Type"] = "PSTH"

rawAccData = t(allFeatures)
colnames(rawAccData) = currAccWnn$uF$cellClass
rownames(rawAccData) = seq(1,nreps)
rawAccDataFeatures = melt(rawAccData)
rawAccDataFeatures["Type"] = "Features"

combData = rbind(rawAccDataWF, rawAccDataWNN, rawAccDataFeatures, rawAccDataPSTH, rawAccDataISI)
colnames(combData) = c("Run","CellType","Acc","Modality")

library(ggthemes)
summaryData = dataSummary(combData, varname="Acc", groupnames = c("CellType","Modality"))
summaryData$se = summaryData$sd/sqrt(nreps)

#summaryData <- rbind(summaryData,averageClassConcat)

p<- ggplot(summaryData, aes(x=CellType, y=Acc, group=Modality, color=Modality))
p = p + theme_classic() + 
  coord_cartesian(clip="off") +
  geom_point(aes(size=4), position=position_dodge(1)) +
  geom_line(aes(size=0.2), position=position_dodge(1)) +
  geom_errorbar(aes(ymin=Acc-se, ymax=Acc+se), width=.2,
                position=position_dodge(1)) 
p = p + theme(text=element_text(size=20)) + ylim(50,100)
p = p + ggtitle(paste0("Classifier at Embedding-D = ",as.character(UMAP.components)))
show(p)

ggsave(paste0(as.character(UMAP.components),".jpg"), width = 10, height = 7)

# p<- ggplot(summaryData, aes(x=CellType, y=Acc))
# p = p + theme_classic() +
#   coord_cartesian(clip="off") +
#   geom_point(aes(size=6), position=position_dodge(1)) +
#   geom_line(aes(size=0.2), position=position_dodge(1)) +
#   geom_errorbar(aes(ymin=Acc-se, ymax=Acc+se), width=.2,
#                 position=position_dodge(1))
# p = p + theme(text=element_text(size=20))

# se <- apply(allAcc, 1, sd)/sqrt(20)
# p2 = p + geom_point(data=test,aes(x=c(1,2,3,4,5),y=concatAcc,size=20)) +
#  ylim(25,100) +
#  geom_errorbar(data=test,aes(x=c(1,2,3,4,5),y=concatAcc,ymin=concatAcc-se,ymax=concatAcc+se))


#weights = data.frame(juxtaData$WF.weight, juxtaData$ISI.weight, juxtaData$PSTH.weight);
#colnames(weights) = c('WF','ISI', "PSTH")
#colMeans(weights)

