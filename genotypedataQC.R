library(tidyverse)
library(dplyr)

geno.qc <- genotypes


names(geno.qc) <- gsub("[A-Za-z-]", "", names(geno.qc)) #remove all letters and dashes from column names 

groups <- split(seq_len(ncol(geno.qc)), colnames(geno.qc)) #which columns have the same name 
groups <- groups[lengths(groups) > 1]

pct_diff <- sapply(groups, function(idx) {
  sub <- geno.qc[, idx, drop = FALSE]
  agree <- apply(sub, 1, function(r) length(unique(r)) == 1) #what is the percent of cels that are differeent between groups? 
  100 * mean(!agree)
})

pct_diff <-round(sort(pct_diff, decreasing = TRUE), 2) #sort 

percent_different <- as.data.frame(pct_diff) #as data frame 


length(percent_different$pct_diff)
sum(percent_different$pct_diff > 10, na.rm = TRUE)/252. #30% are more than 10% different 

