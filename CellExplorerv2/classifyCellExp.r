
doClassify = function(E, data, numIter = 10, numreps=2, 
                      method='boot')
{
  tempCells = str_trim(data$cType)
  idx = tempCells %in% c("Pyra", "PV", "SST","Juxt","Axo_","VIP")
  # idx = tempCells %in% c("Pyra","PV", "SST")
  

  E = E[idx,]
  dim(E)
  tempCells = tempCells[idx]
  origCells = factor(tempCells);
  E$origCells = origCells
  
  Rs = runif(numIter);
  D = c()
  AccO = c()
  for(k in 1:numIter)
  {
   set.seed(k)
    print(k)
    i <- createDataPartition(E$origCells, times = 1, p = 0.7, list = FALSE)
    training = E[i[,1],]
    testingset = E[-i[,1],]
   
    ctrl <- trainControl(method = "repeatedcv", number=numreps)
    #fit a regression model and use k-fold CV to evaluate performance
    model <- train(origCells~., data = training, method = "rbf", 
                   trControl = ctrl, verbose=FALSE)
    mean(model$results$Accuracy)
    cVpc = predict(model, newdata = testingset)
    cM = confusionMatrix(predict(model, newdata = testingset),testingset$origCells)
    X = data.frame(cM$byClass)
    X$Balanced.Accuracy
    
    #X = X[c(4, 2, 3, 5, 1, 6),]
    X = X[c(4, 3, 5,1,2, 6),]
    
    D = rbind(D, X$Balanced.Accuracy)
    AccO[k] = cM$overall[1]
  }


  uF = data.frame(AccV = D*100)
  colnames(uF) = rownames(X)
  return(list(uF = uF, AccO = AccO))
}

postProcess = function(accData, type)
{
  currAcc = accData$uF
  nreps = dim(currAcc)[1]
  rawAccData = currAcc
  rownames(currAcc) = seq(1,nreps)
  currAccOut = melt(currAcc)
  currAccOut["Type"] = type
  return(currAccOut)
}


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

data <- RunUMAP(data, nn.name = "weighted.nn", 
                reduction.name = "wnn.umap2", 
                reduction.key = "wnnUMAP2_", seed.use=UMAP.SEED, 
                n.components = 10
)




E1 = data.frame(Embeddings(data[["wnn.umap2"]]))
E2 = data.frame(Embeddings(data[["WFumap"]]))
E3 = data.frame(Embeddings(data[["ISI1umap"]]))
E4 = data.frame(Embeddings(data[["featuresumap"]]))
E5 = data.frame(Embeddings(data[["ACGumap"]]))




AccWnn = doClassify(E1, data)
rawAccWNN = postProcess(AccWnn, "WNN")

AccWf = doClassify(E2, data)
rawAccWf = postProcess(AccWf, "WF")

AccISI = doClassify(E3, data)
rawAccISI = postProcess(AccISI, "ISI1")

AccFeatures = doClassify(E4, data)
rawAccfeatures = postProcess(AccFeatures, "features")

AccACG = doClassify(E5, data)
rawAccACG = postProcess(AccACG, "ACG")

combData = rbind(rawAccWNN, rawAccWf, rawAccISI, rawAccfeatures, rawAccACG)

AvgData = data.frame(Wnn = AccWnn$AccO, Wf = AccWf$AccO, ACG = AccACG$AccO, isi = AccISI$AccO, features
                     =AccFeatures$AccO)



colnames(combData) = c("CellType","Acc","Modality")
nreps = 20
library(ggthemes)
summaryData = dataSummary(combData, varname="Acc", groupnames = c("CellType","Modality"))
summaryData$se = summaryData$sd/sqrt(nreps)
p<- ggplot(summaryData, aes(x=CellType, y=Acc, group=Modality, color=Modality)) + 
  geom_line() + theme_minimal() + 
  geom_point(aes(size=8)) +
  geom_errorbar(aes(ymin=Acc-se, ymax=Acc+se), width=.2,
                position=position_dodge(0.0))

show(p)
