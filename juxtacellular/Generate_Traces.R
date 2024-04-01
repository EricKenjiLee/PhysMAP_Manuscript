rm(list=ls()) 

library(reshape)
library(ggplot2)

basedir <- dirname(sys.frame(1)$ofile)
setwd(basedir)

source("./helperFunctions.r")

allData = readMat("./JianingData/MergedData.mat");

plot.waveforms = function(data,layer.cell.type,modality)
{ layerData = data$layer;
  cellType = data$true.cell.type;
  layerCellType = paste0(trimws(cellType),'-', trimws(layerData));

  if (modality == "waveform"){
    unit_data <- data$waveformV[grepl(paste(layer.cell.type,collapse="|"),layerCellType),]
  } else if (modality == 'ISI'){
    unit_data <- data$ISI[grep(layer.cell.type,layerCellType),]
  } else if (modality == 'PSTH'){
    unit_data <- data$AllResp[grep(layer.cell.type,layerCellType),]
  }
  
  
  colnames(layerData) = "Layer";
  
  colnames(cellType) = "CellType"
  
  units = data.frame(unit_data)
  avg_trace_df <- data.frame(colMeans(unit_data))
  units <- melt(units)
  units$rowid <- 1:nrow(unit_data)
  # 
  # mean_data <- units %>%
  #   group_by(factor(rowid)) %>%
  #   summarise(mean_value = mean(value))
  # 
  p <- ggplot(data = units) +
    geom_line(data=units,aes(x = variable, y = value, group=factor(rowid), color = factor(rowid)),size=1, alpha=0.3, show.legend=FALSE) + theme_void() + 
    scale_color_manual(values=c(rep('gray',nrow(units)))) + 
    geom_line(data = avg_trace_df, aes(x = 1:length(colMeans.unit_data.),y=colMeans.unit_data.),size=2)
  return(p)
}


