setwd("~/Documents/usda potato breeding/data_analysis/raw_data/gravity2023")
#column names to vector 
col_names <- c("zero", "Sample", "DryWeightTotal", "WetWeigthTotal", "DryWeight", "WetWeight", "SpecificGravity", "Date", "Time", "User")
getwd()

#write file names to vector
files <- list.files(pattern = "\\.csv$")
#empty container 
all_data <- data.frame()


for (f in files) {
  #read file into dataframe, apply column name from vector
  df <- read.csv(f, header = FALSE, col.names = col_names)
  
  df$date <- as.Date(substr(f, 1, 10))   # apply file name to date column as date (file name should be date YYYY-MM-DD)
  all_data <- rbind(all_data, df)
}

write.csv(all_data, "combined_data.csv", row.names = FALSE)

