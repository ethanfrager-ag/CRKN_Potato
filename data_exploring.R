library(tidyverse)
library(ggplot2)
#each section is an independent step ---

#read in csv
combined <- read.csv("~/Documents/usda potato breeding/data_analysis/clean_data/gravity2324.csv", stringsAsFactors = FALSE)     
#remove "zero" column 
combined$zero <- NULL
write.csv(combined, "combined_2024.csv", row.names = FALSE)

#find high sg values (impossible ones)
sum(combined$NewSG > 1.2)
combined$NewSG[combined$NewSG > 1.2]
which(combined$NewSG > 1.2)
#find low SG
sum(combined$NewSG < 0.8)
combined$NewSG[combined$NewSG < 0.8]
which(combined$NewSG < 0.8)

#find how many entries have a SG in ideal range
sgmatrix <- as.matrix(combined$NewSG)
idealSGs <- sgmatrix[sgmatrix >= 1.075 & sgmatrix <= 1.08]
length(idealSGs)

#ecdf plot of sg 
plot(ecdf(combined$NewSG), xlim = c(0.85, 1.25), main = "cumulative proportion of SG values", xlab = "Specific Gravity")
abline(h = c(0.25, 0.5, 0.75), col = "grey70", lty = 3)

#sg density plot with rug 
plot(density(combined$NewSG, na.rm = TRUE), xlim = c(0.85, 1.25), main = "")
rug(combined$NewSG, col = rgb(0, 0, 0, 0.3))

#combine 2024 and 2023 csv

setwd("~/Documents/usda potato breeding/data_analysis/clean_data")

a <- read.csv("combined_2023.csv", stringsAsFactors = FALSE)
b <- read.csv("combined_2024.csv", stringsAsFactors = FALSE)

setdiff(names(a), names(b))   # in a but not b
setdiff(names(b), names(a))   # in b but not a
combined_clean <- rbind(a, b)
write.csv(combined_clean, "gravity2324.csv", row.names = FALSE)

#replace NA with CONTROL 
df <- read.csv("gravity2324.csv")
  
df[is.na(df)] <- "CONTROL"


write.csv(df, "gravity2324.csv", row.names = FALSE)
