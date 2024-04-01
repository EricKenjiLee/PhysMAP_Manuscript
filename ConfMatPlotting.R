setwd("/Users/kenjilee/Downloads/conf_mats")

# Get a list of all .Rda files in the directory
file_list <- list.files(pattern = "\\.Rda$", full.names = TRUE)

# Create an empty list to store the loaded data frames
data_list <- list()

# Iterate through each file and load the data into the list
i = 1
for (file in file_list) {
  # Load the data from the .Rda file
  load(file)
  loaded_data <- conf_mat$`Confusion Matrix`[[1]][,"N"]

  
  # Assuming the data frames have the same structure, you can perform any additional processing here
  
  # Add the loaded data to the list
  data_list[[i]] <- loaded_data
  i = i + 1
}

confMatAvg = data_list[[1]]
for (i in 2:length(data_list)){
  confMatAvg = confMatAvg + data_list[[i]]
}

confMatAvg = confMatAvg/10

avg_conf_mat$`Confusion Matrix`[[1]] = conf_mat$`Confusion Matrix`[[1]]
avg_conf_mat$`Confusion Matrix`[[1]]["N"] = confMatAvg

p_conf <- plot_confusion_matrix(
  avg_conf_mat$`Confusion Matrix`[[1]],
  intensity_by='normalized',
  palette="Blues",add_counts=FALSE,add_normalized=FALSE,
  add_arrows = FALSE
)
