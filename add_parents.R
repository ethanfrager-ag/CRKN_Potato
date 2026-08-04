#open files
combined <- read.csv("~/Documents/usda potato breeding/data_analysis/raw_data/gravity2023/combined_data.csv", stringsAsFactors = FALSE)     
Field    <- read.csv("~/Documents/usda potato breeding/data_analysis/raw_data/metadata2023/CRKNOthelloFieldKey_2023.csv", stringsAsFactors = FALSE)            
FamilyN  <- read.csv("~/Documents/usda potato breeding/data_analysis/raw_data/metadata2023/CRKNFamilyNameKey.csv", stringsAsFactors = FALSE)            
#run first 
# last 4 characters of combined$Sample, as an integer
s <- as.character(combined$Sample)
key <- as.integer(substr(s, nchar(s) - 3, nchar(s)))

# pull the matching Clone code out of Field
combined$Clone <- Field$Clone[match(key, Field$Sample)]

write.csv(combined, "combined.csv", row.names = FALSE)

#Step 2 (run after)

combined <- read.csv("~/Documents/usda potato breeding/data_analysis/raw_data/gravity2023/combined.csv", stringsAsFactors = FALSE)     
#remove hyphen and beyond
combined$FamilyN <- sub("-.*$", "", as.character(combined$Clone))

ss <- as.character(combined$FamilyN)
key <- as.integer(substr(ss, nchar(ss) - 3, nchar(ss)))
combined$Male <- FamilyN$Male[match(combined$FamilyN, FamilyN$Field.name)]
combined$Female <- FamilyN$Female[match(combined$FamilyN, FamilyN$Field.name)]


write.csv(combined, "combined.csv", row.names = FALSE)
