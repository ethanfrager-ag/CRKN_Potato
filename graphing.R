library(car)
library(ggplot2)

#qq plot 
qqnorm(weightsag$Gravity, main = "Specific Gravity"); qqline(weightsag$Gravity, col = "red")
shapiro.test(combined$AvgSG)
qqnorm(sqrt(weightsag$DryWeight), main = "SQRT Dry Weight of Tubers"); qqline(sqrt(weightsag$DryWeight), col = "red")

#single t-test
t.test(weightsag$Gravity, mu = 1.0751)
t.test(weightsag$DryWeight, mu = 8.54)

#new vectors SG and Yield by year 
sg2025 <- combined$AvgSG[combined$Year == 2025] 
sg2024 <- combined$AvgSG[combined$Year == 2024] 
sg2023 <- combined$AvgSG[combined$Year == 2023] 

dryw2025 <- combined$DryAvg[combined$Year == 2025] #dry weight 
dryw2024 <- combined$DryAvg[combined$Year == 2024] #dry weight 
dryw2023 <- combined$DryAvg[combined$Year == 2023] #dry weight 




#table of dry weight by year 
cnt <- aggregate(DryWeight ~ Sample + Year, data = combined, FUN = length)
table(cnt$Year, cnt$DryWeight)


#boxplot of sg x year + p values

library(tidyverse)


P <- weightsag %>%
    ggplot( aes(Year, Gravity, Year)) +
    geom_boxplot() +
  
    ggtitle("Dry Weight By Year") +
    xlab("Year")

P + stat_compare_means()

compare_means(Gravity ~ Year,  data = weightsag)

my_comparisons <- list( c("2023", "2024") )
ggboxplot(weightsag, x = "Year", y = "Gravity",
          color = "Year", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1)     # Add global p-value


#venn diagram of how many clones overlap by year 

library(ggVennDiagram)

sg2023 <- unique(weightsag$Clone[weightsag$Year == "2023"]) #new vectors by year 
sg2024 <- unique(weightsag$Clone[weightsag$Year == "2024"])

sg2023 <- sg2023[!grepl("CONTROL", sg2023)] #remove "CONTROL"
sg2024 <- sg2024[!grepl("CONTROL", sg2024)]

clone_list <- list(sg2023, sg2024)


ggVennDiagram(clone_list) +
  scale_fill_gradient(low = "grey90", high = "#534AB7") +
labs(title = "Clone overlap by year")


sum(table(unique(weightsag$Clone[weightsag$Year == "2023"]))) #double check
sum(table(unique(weightsag$Clone[weightsag$Year == "2024"])))
#plot of how many clones per year 
yearcount <- table(weightsag$Year)

barplot(yearcount) + title("Clones Per Year ")

#plot dry weight x sg 

cor.test(weightsag$DryWeight, weightsag$Gravity)
plot(weightsag$DryWeight, weightsag$Gravity,
     xlab = "dry weight", ylab = "specific gravity", main = "cor = 0.1849709")

#yield by user


usercount <- table(weightsag$User) #measurements per user

barplot(usercount) + title("Measurements per User ") 

  userlist <- list(
    "AP" = unique(weightsag$DryWeight[weightsag$User == "AP"]),
    "IB" = unique(weightsag$DryWeight[weightsag$User == "IB"]),
    "NM" = unique(weightsag$DryWeight[weightsag$User == "NM"]),
    "CH" = unique(weightsag$DryWeight[weightsag$User == "CH"]),
    "CR" = unique(weightsag$DryWeight[weightsag$User == "CR"]),
    "KR" = unique(weightsag$DryWeight[weightsag$User == "KR"]),
    "VP" = unique(weightsag$DryWeight[weightsag$User == "VP"])

  )
boxplot(userlist, xlab = "User", ylab = "Dry weight",
        main = "Dry weight by User")

#sg by user box plot

userlistG <- list(
  "AP" = unique(weightsag$Gravity[weightsag$User == "AP"]),
  "IB" = unique(weightsag$Gravity[weightsag$User == "IB"]),
  "NM" = unique(weightsag$Gravity[weightsag$User == "NM"]),
  "CH" = unique(weightsag$Gravity[weightsag$User == "CH"]),
  "CR" = unique(weightsag$Gravity[weightsag$User == "CR"]),
  "KR" = unique(weightsag$Gravity[weightsag$User == "KR"]),
  "VP" = unique(weightsag$Gravity[weightsag$User == "VP"])
  
)
boxplot(userlistG, xlab = "User", ylab = "Specific Gravity",
        main = "Specific Gravity by User")

# plot to correlate Sg 2024 x SG 2023 and find Rsquared

d23 <- weightsag[weightsag$Year == "2023", c("Clone", "Gravity")]
d24 <- weightsag[weightsag$Year == "2024", c("Clone", "Gravity")]
names(d23)[2] <- "SG_2023"
names(d24)[2] <- "SG_2024"

paired <- merge(d23, d24, by = "Clone")

plot(paired$SG_2023, paired$SG_2024, xlab = "SG 2023", ylab = "SG 2024")
abline(0, 1, col = "grey", lty = 2)   # reference line: equal in both years
legend("topleft", legend = "R^2 = 0.300", bty = "n")

fit <- lm(SG_2024 ~ SG_2023, data = paired)
abline(fit, col = "red")
summary(fit)$r.squared
# plot to correlate dw 2024 x dw 2023 and find Rsquared

d23 <- weightsag[weightsag$Year == "2023", c("Clone", "DryWeight")]
d24 <- weightsag[weightsag$Year == "2024", c("Clone", "DryWeight")]
names(d23)[2] <- "SG_2023"
names(d24)[2] <- "SG_2024"

paired <- merge(d23, d24, by = "Clone")

plot(paired$SG_2023, paired$SG_2024, xlab = "Yield 2023", ylab = "Yield 2024")
abline(0, 1, col = "grey", lty = 2)   # reference line: equal in both years
legend("topleft", legend = "R^2 = 0.123", bty = "n")

fit <- lm(SG_2024 ~ SG_2023, data = paired)
abline(fit, col = "red")
summary(fit)$r.squared

#histogram of rank change yield
hist(comp$RankChange, xlab = "Rank Change", main = "Change in Yield Ranks Between 2023 and 2024", col = "grey")
qqnorm(comp$RankChange)

#histogram of rank change sg
hist(compsg$RankChange, xlab = "Rank Change", main = "Change in SG Ranks Between 2023 and 2024", col = "grey")
qqnorm(compsg$RankChange)

#plot anova of yield

library(ggstatsplot)

ggbetweenstats(
  data = onlyclonesanova,
  x = FamilyN,
  y = DryWeight,
  xlab = "Control Variety",
  ylab = "Yield in Kg",
  type = "parametric", # ANOVA or Kruskal-Wallis
  var.equal = FALSE, # ANOVA or Welch ANOVA
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)

