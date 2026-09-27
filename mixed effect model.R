library(emmeans)
library(tidyverse)
library(car)
library(dplyr)
library(multcomp)      
library(multcompView)


###### Creating a fixed effect linear model of how year, male, female, and male + female affect yield
#remove CONTROL entries
#wt <- droplevels(subset(weightsag, Male != "CONTROL" & Female != "CONTROL"))
wt <- weightsag

# 9/23/26: update "wt" with pedigree corrected using trios 
wt.new <- wt
clean <- function(x) trimws(as.character(x))

ped <- triopedigree[!duplicated(clean(triopedigree$id)), ]          # match() takes the first hit only
i   <- match(clean(wt.new$Clone), clean(ped$id))      # NA where the clone isn't listed

new_sire <- clean(ped$father)[i]
new_dam  <- clean(ped$mother)[i]

ok_s <- !is.na(new_sire) & new_sire != ""         # don't overwrite with blanks
ok_d <- !is.na(new_dam)  & new_dam  != ""

wt.new$Male[ok_s] <- new_sire[ok_s]
wt.new$Female[ok_d] <- new_dam[ok_d]

cat("Matched:", sum(!is.na(i)), "of", nrow(wt.new), "rows\n")
######


levels(wt$Male)[levels(wt$Male) == "POR08BD3-1"] <- "POR08BD1-3" #correct wrong name 
#remove control and male outlier 
wt <- subset(wt, Male != "AO02183-2" &
               Clone != "CONTROL" &
               Male  != "CONTROL" &
               Female != "CONTROL")
wt <- droplevels(wt)
wt$Male[wt$Male == "POR08BD3-3"] <- "POR08BD1-3"
wt.new$Male[wt$Male == "POR08BD3-1"] <- "POR08BD1-3"



modelyield <- lm(sqrt(DryWeight) ~ Female*Male + Year, data = wt.new) #model
summary(modelyield)
qqnorm(resid(modelyield))

#graph anova
modelyieldanova <- as.data.frame(Anova(modelyield, type = 2)) 
modelyieldanova$Term <- rownames(modelyieldanova)
modelyieldanova$PctVar <- modelyieldanova$`Sum Sq` / sum(modelyieldanova$`Sum Sq`) * 100

ggplot(modelyieldanova, aes(x = reorder(Term, PctVar), y = PctVar)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "% of variance in yield explained")
#compute eta squared
eta_squared(Anova(modelyield, type = 2), partial = TRUE)
#descrptive statustics for anova table 
wt %>%
  group_by(Male) %>%
  summarise(
    n      = n(),
    mean   = mean(DryWeight, na.rm = TRUE),
    sd     = sd(DryWeight, na.rm = TRUE),
    .groups = "drop"
  )

#plot estimated marginal means linear model
em <- as.data.frame(emmeans(modelyield, ~ Male, type = "response"))
ggplot(na.omit(em), aes(x = reorder(Male, response), y = response)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, na.rm = TRUE) +
  coord_flip() +
  labs(x = "Male parent", y = "Estimated marginal mean (Yield in Kg)")
###### linear model of total starch in grams
#square root transform to make the data more normal. 
modelsolids <- lm(sqrt(TotalStarchGrams) ~ Year + Male*Female, data = wt.new)
#modelsolids <- lmer(solids ~ Year + Male + Female + (1|Clone), data = wt)

summary(modelsolids)
qqnorm(resid(modelsolids))



em <- as.data.frame(emmeans(modelsolids, ~ Female))
ggplot(na.omit(em), aes(x = reorder(Female, emmean), y = emmean,)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, na.rm = TRUE) +
  coord_flip() +
  labs(x = "Female Parent", y = "Estimated Marginal Means (Starch In Grams)")


#graph anova of solids emmeans
