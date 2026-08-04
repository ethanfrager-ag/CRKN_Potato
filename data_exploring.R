library(tidyverse)
#read in csv
combined <- read.csv("~/Documents/usda potato breeding/data_analysis/raw_data/gravity2023/combined.csv", stringsAsFactors = FALSE)     
#remove "zero" column 
combined$zero <- NULL
write.csv(combined, "combined.csv", row.names = FALSE)

#find high sg values (impossible ones)
sum(combined$SpecificGravity > 1.2)
combined$SpecificGravity[combined$SpecificGravity > 1.2]
which(combined$SpecificGravity > 1.2)
#find low SG
sum(combined$SpecificGravity < 0.8)
combined$SpecificGravity[combined$SpecificGravity < 0.8]
which(combined$SpecificGravity < 0.8)

#find how many entries have a SG in ideal range
sgmatrix <- as.matrix(combined$SpecificGravity)
idealSGs <- sgmatrix[sgmatrix >= 1.075 & sgmatrix <= 1.08]
length(idealSGs)