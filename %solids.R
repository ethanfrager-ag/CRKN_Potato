weightsag <- read.csv("~/Documents/usda potato breeding/data_analysis/clean_data/combined_2024.csv", stringsAsFactors = FALSE)     


b = -214.9206
m = 218.1852

weightsag$solids <- (m*weightsag$Gravity)+b

write.csv(weightsag, "agdata.csv", row.names = FALSE)

#equation from: Kleinkopf, G. E., et al. "Specific gravity of Russet Burbank potatoes." American potato journal 64.11 (1987): 579-587.

#recalculate SG from wet and dry weight 
setwd("~/Documents/usda potato breeding/data_analysis/clean_data")

weightsag$Gravity <- (weightsag$DryWeight/(weightsag$DryWeight-weightsag$WetWeight))

write.csv(weightsag, "agdata.csv", row.names = FALSE)

#find total weight of starch 

weightsag$TotalStarchGrams <- (((0.01*weightsag$solids)*weightsag$DryWeight)*1000)
write.csv(weightsag, "agdata.csv", row.names = FALSE)


