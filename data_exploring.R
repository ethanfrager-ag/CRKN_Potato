library(tidyverse)
library(ggplot2)
#each section is an independent step ---

#read in csv
combined <- read.csv("~/Documents/usda potato breeding/data_analysis/clean_data/gravity2324.csv", stringsAsFactors = FALSE)  
weightsag <- read.csv("~/Documents/usda potato breeding/data_analysis/R/agdata.csv", stringsAsFactors = FALSE)     

#remove "zero" column 
combined$zero <- NULL
write.csv(combined, "combined_2024.csv", row.names = FALSE)

#find high sg values (impossible ones)
sum(weightsag$Gravity > 1.125)
weightsag$Gravity[weightsag$Gravity > 1.125]
which(weightsag$Gravity > 1.125)
weightsag <- weightsag[-which(weightsag$Gravity > 1.125), ] #remove values > 1.125

#find low SG
sum(weightsag$Gravity < 1.02)
weightsag$Gravity[weightsag$Gravity < 1.02]
which(weightsag$Gravity < 1.02)
weightsag <- weightsag[-which(weightsag$Gravity < 1.02), ] #remove values < 1.02

#find how many entries have a SG in ideal range
sgmatrix <- as.matrix(weightsag$Gravity)
idealSGs <- sgmatrix[sgmatrix >= 1.075 & sgmatrix <= 1.08]
length(idealSGs)

#ecdf plot of sg 
plot(ecdf(combined$NewSG), xlim = c(0.85, 1.25), main = "cumulative proportion of SG values", xlab = "Specific Gravity")
abline(h = c(0.25, 0.5, 0.75), col = "grey70", lty = 3)

#sg density plot with rug 
plot(density(combined$AvgSG, na.rm = TRUE), xlim = c(0.85, 1.25), main = "")
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



#aggregate all data into new csv with no reapeated sample codes


weightsag <- aggregate(cbind(DryWeight, WetWeight) ~ Sample + Year,
                       data = combined, FUN = sum)

info <- combined[!duplicated(combined$Sample),
                 c("Sample", "File.date", "Date", "User", "Time",
                   "Clone", "FamilyN", "Male", "Female")]

weightsag <- merge(weightsag, info, by = "Sample", sort = FALSE)
weightsag <- weightsag[order(weightsag$Year), ]
write.csv(weightsag,"agdata.csv", row.names = FALSE)

#rank starch content in grams and put in new df
starchrank <- weightsag[order(-weightsag$TotalStarchGrams), ]


starchtop10 <- head(starchrank, n = 10)
write.csv(starchtop10, "starchtop10.csv", row.names = FALSE)



plot(weightsag$Gravity, weightsag$TotalStarchLBs, type = "l")
hist(starchrank$TotalStarchGrams)

#correct error in year column 
weightsag$Year <- as.character(weightsag$Year) #aggregated data
weightsag$Year <- ifelse(weightsag$Year == "2025", 2024, 2023)
weightsag$Year <- factor(weightsag$Year)
write.csv(weightsag, "agdata.csv", row.names = FALSE)

combined$Year <- as.character(combined$Year) #unaggregated data
combined$Year <- ifelse(combined$Year == "2025", 2024, 2023)
combined$Year <- factor(combined$Year)
write.csv(combined, "gravity2324.csv", row.names = FALSE)

# 2 sample t tests 


tapply(weightsag$Gravity, weightsag$Year, summary) #summary of both years 
tapply(weightsag$DryWeight, weightsag$Year, summary)


tapply(weightsag$Gravity, weightsag$Year, sd) #standard deviation in both years 
tapply(weightsag$DryWeight, weightsag$Year, sd)

t.test(Gravity ~ Year, data = weightsag) #ttests
t.test(DryWeight ~ Year, data = weightsag)


#find and compare rank order of Yield in both years
d2024 <- weightsag[weightsag$Year == "2024", ] #vectors by year
d2023 <- weightsag[weightsag$Year == "2023", ]

rank2024 <- d2024[order(-d2024$DryWeight), ] #order vectors by dry weight 
rank2023 <- d2023[order(-d2023$DryWeight), ]

rank2024$Rank2024 <- seq_len(nrow(rank2024)) #make columns with the rank assigned 
rank2023$Rank2023 <- seq_len(nrow(rank2023))

comp <- merge(rank2023[, c("Clone", "Rank2023")], #combine the two vectors and remove duplicates
              rank2024[, c("Clone", "Rank2024")],
              by = "Clone")

comp$RankChange <- comp$Rank2024 - comp$Rank2023 #column with how the rank changes

sum(comp$Rank2023 != comp$Rank2024)          # any change at all
sum(abs(comp$RankChange) > 5)                # moved more than 5 places

cor(comp$Rank2023, comp$Rank2024) #correlation test 
comp[order(-abs(comp$RankChange)), ][1:10, ]
summary(comp$RankChange)

#find and compare rank order of specific gravity  in both years
dg2024 <- weightsag[weightsag$Year == "2024", ] #vectors by year
dg2023 <- weightsag[weightsag$Year == "2023", ]

ranksg2024 <- d2024[order(-d2024$Gravity), ] #order vectors by SG
ranksg2023 <- d2023[order(-d2023$Gravity), ]

ranksg2024$Rank2024 <- seq_len(nrow(ranksg2024)) #now that they are ordered, make a column with n+1 until the end of sequence (this is the rank #)
ranksg2023$Rank2023 <- seq_len(nrow(ranksg2023))

compsg <- merge(ranksg2023[, c("Clone", "Rank2023")], #combine the two vectors and remove duplicates
              ranksg2024[, c("Clone", "Rank2024")],
              by = "Clone")

compsg$RankChange <- compsg$Rank2024 - compsg$Rank2023 #column with how the rank changes

sum(compsg$Rank2023 != compsg$Rank2024)          # any change at all
sum(abs(compsg$RankChange) > 5)                # moved more than 5 places

cor(compsg$Rank2023, compsg$Rank2024) #correlation test 
compsg[order(-abs(compsg$RankChange)), ][1:10, ] #top 10 biggest changes in rank
summary(compsg$RankChange)


#anova to propose cut off value for yield
library(car)
onlyclonesanova <- weightsag
onlyclonesanova <- sub("_.*$", "", as.character(onlyclonesanova$FamilyN)) #remove everything in FamilyN past the _

targetclones <- c("Burbank", "Ivory", "Alturus","Umatilla","Ranger") #List of check varieties 
onlyclonesanova <- weightsag[weightsag$FamilyN %in% targetclones, ] #move only target clones to df 

factor(onlyclonesanova$FamilyN) #factorize family name 
fit <- aov(DryWeight ~ FamilyN, data = onlyclonesanova, var.equal = FALSE) #welch anova type 1 
pairwise.t.test(onlyclonesanova$DryWeight, onlyclonesanova$FamilyN,   #must use holm pairwise t-test instead of tukey due to welch anova
                p.adjust.method = "holm"
)
TukeyHSD(fit)
print(fit)

qqnorm(fit$residuals); qqline(fit$residuals, col = "red") #testing 4 normality of residuals
leveneTest(DryWeight ~ FamilyN, data = onlyclonesanova, var.equal = FALSE) #homogeneity of varriance? (not very conclusive)
boxplot(DryWeight ~ FamilyN, data = onlyclonesanova)
tapply(onlyclonesanova$DryWeight, onlyclonesanova$FamilyN, sd) #find sd & how many entries
tapply(onlyclonesanova$DryWeight, onlyclonesanova$FamilyN, length)

aggregate(DryWeight ~ FamilyN,  #descriptive stats
          data = onlyclonesanova,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)


#plot anova of yield

library(ggstatsplot)

ggbetweenstats(
  data = onlyclonesanova,
  x = FamilyN,
  y = DryWeight,
  xlab = "Control Variety",
  ylab = "Yield in Kg",
  type = "parametric", # ANOVA or Kruskal-Wallis
  var.equal = TRUE, # ANOVA or Welch ANOVA
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)



#anova to propose cut off value for SG
fit_sg <- aov(Gravity ~ FamilyN, data = onlyclonesanova) # anova type 1 
TukeyHSD(fit_sg)
print(fit_sg)

qqnorm(fit_sg$residuals); qqline(fit_sg$residuals, col = "red") #testing 4 normality of residuals
leveneTest(Gravity ~ FamilyN, data = onlyclonesanova, var.equal = FALSE) #homogeneity of varriance?
boxplot(Gravity ~ FamilyN, data = onlyclonesanova)
tapply(onlyclonesanova$Gravity, onlyclonesanova$FamilyN, sd) #sd is not as variable as with yield so regular anova
tapply(onlyclonesanova$Gravity, onlyclonesanova$FamilyN, length)

aggregate(Gravity ~ FamilyN,  #descriptive stats
          data = onlyclonesanova,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)


#plot anova of sg

library(ggstatsplot)

ggbetweenstats(
  data = onlyclonesanova,
  x = FamilyN,
  y = Gravity,
  xlab = "Control Variety",
  ylab = "Specific Gravity of Tubers",
  type = "parametric", # ANOVA or Kruskal-Wallis
  var.equal = TRUE, # ANOVA or Welch ANOVA
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)





