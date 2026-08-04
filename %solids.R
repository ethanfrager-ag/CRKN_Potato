data2023 <- read.csv("~/Documents/usda potato breeding/data_analysis/raw_data/gravity2023/combined.csv", stringsAsFactors = FALSE)     


b = -214.9206
m = 218.1852

data2023$solids <- (m*data2023$SpecificGravity)+b

write.csv(data2023, "combined.csv", row.names = FALSE)

#equation from: Kleinkopf, G. E., et al. "Specific gravity of Russet Burbank potatoes." American potato journal 64.11 (1987): 579-587.
