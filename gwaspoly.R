library(GWASpoly)
VCF2dosage("~/Documents/usda potato breeding/data_analysis/R/DArTag.vcf.gz","./gwasdosage.csv", "GT", 4, samples=NULL,
           min.DP=1, 0.1, min.minor=5) #convert vcf to dosage file 
gwasdosage <- read.csv("gwasdosage.csv") #read in the dosage file 
#format the dosagefile data
colnames(gwasdosage) <- sub("^[^_]*_([^_]+)_(\\d+).*$", "\\1-\\2", colnames(gwasdosage)) #change the format to same as phenofile 
#remove 3 malformed rows
gwasdosage <- gwasdosage[-c(3855,3856,3857),]
write.csv(gwasdosage, "gwasdosage.csv", row.names = FALSE)



phenofilegwas <- data2324[c("Clone","DryWeight","SpecificGravity","Year")]
names(phenofilegwas)[4] <- "env"
names(phenofilegwas)[1] <- "id"
phenofilegwas$env <- factor(phenofilegwas$env)



write.csv(phenofilegwas, "phenofilegwas.csv", row.names = FALSE)
#make the phenotype file
#year is fixed effect, traits are yeild and sg 

gwasdata <- read.GWASpoly(ploidy=4, pheno.file="phenofilegwas.csv", geno.file="gwasdosage.csv",
                          format="numeric", n.traits=2, delim=",") #read the data into gwas poly (dosage and phenotype) 

data.LOCO <- set.K(gwasdata, n.core = 2, LOCO = TRUE) #apply the K model without LOCO to control for population structure 


Npop <- 843 #Population size
params <- set.params(geno.freq = 1 - 5/Npop, fixed = "env", fixed.type = "factor") #parameters (remove markers below frequency of 0.05)
#compute additive and single domoinance p values for all markers (basically the whole gwas)
data.loco.scan <- GWASpoly(data=data.LOCO,models=c("additive","1-dom"),
                           traits=c("DryWeight","SpecificGravity"),params=params,n.core=2) 


#it worked
#qq plotr 
library(ggplot2)
qq.plot(datagwas2,trait="DryWeight") + ggtitle(label="LOCO yield")
qq.plot(datagwas2,trait="SpecificGravity") + ggtitle(label="LOCO gravity")
#set a threshhold of 0.05 to control fasle positive rate 
datagwas2 <- set.threshold(data.loco.scan,method="M.eff",level=0.05) 

#there is an error because of monomorphic markers 

#manhatten plots
p <- manhattan.plot(datagwas2,traits="DryWeight")
p + theme(axis.text.x = element_text(angle=90,vjust=0.5))


p <- manhattan.plot(datagwas2,traits="SpecificGravity")
p + theme(axis.text.x = element_text(angle=90,vjust=0.5))

#zoom in on chromosome 5 since it seems to have big hits for yield and sg 
manhattan.plot(datagwas2,traits="SpecificGravity",chrom="chr05")
manhattan.plot(datagwas2,traits="DryWeight",chrom="chr08")
#they seem to hit on exactly the same region 




#find the QTLs
p <- LD.plot(datagwas2, max.loci=1000)
p + xlim(0,100) 

library(knitr)
qtl <- get.QTL(data=datagwas2,traits=c("DryWeight","SpecificGravity"),models="1-dom",bp.window=5e6)
knitr::kable(qtl)

fit.ans <- fit.QTL(data=datagwas2,trait="DryWeight",
                   qtl=qtl[,c("Marker","Model")],
                   fixed=data.frame(Effect="env",Type="factor"))
knitr::kable(fit.ans,digits=3)


#let's run the GWAS on each year individually with no fixed effect. 

#you can go through and change the names to reflect which year you're running on 


phenofilegwas24 <- data2324[c("Clone","DryWeight","SpecificGravity","Year")]
names(phenofilegwas24)[1] <- "id"
phenofile2024 <- phenofilegwas24[phenofilegwas24$Year != 2023, ]
phenofile2024$Year <- NULL

write.csv(phenofile2024, "phenofilegwas2024.csv", row.names = FALSE)

gwasdata24 <- read.GWASpoly(ploidy=4, pheno.file="phenofilegwas2024.csv", geno.file="gwasdosage.csv",
                          format="numeric", n.traits=2, delim=",") #read the data into gwas poly (dosage and phenotype) 

data.LOCO24 <- set.K(gwasdata24, n.core = 2, LOCO = TRUE) #apply the K model with LOCO to control for population structure 


Npop <- 817 #Population size based on the read.GWAS.poly function
params <- set.params(geno.freq = 1 - 5/Npop) #parameters (remove markers below frequency of 0.05). There's no fixed effect. MINOR allele frequency minimum is 0.05

#params <- set.params(geno.freq = 1 - 5/Npop, MAF = 0.05, n.PC = 3) #parameters (remove markers below frequency of 0.05). There's no fixed effect. MINOR allele frequency minimum is 0.05
#compute additive and single domoinance p values for all markers (basically the whole gwas)
data.loco.scan <- GWASpoly(data=data.LOCO24,models=c("additive","1-dom"),
                           traits=c("DryWeight","SpecificGravity"),params=params,n.core=2) 

#grqaphing section

#set a threshhold of 0.05 to control fasle positive rate 
datagwas2 <- set.threshold(data.loco.scan,method="M.eff",level=0.05) 

library(ggplot2)
qq.plot(datagwas2,trait="DryWeight") + ggtitle(label="LOCO yield 2024")
qq.plot(datagwas2,trait="SpecificGravity") + ggtitle(label="LOCO gravity 2024")


#there is an error because of monomorphic markers 

#manhatten plots
p <- manhattan.plot(datagwas2,traits="DryWeight")
p + theme(axis.text.x = element_text(angle=90,vjust=0.5)) + ggtitle(label=" 2024")


p <- manhattan.plot(datagwas2,traits="SpecificGravity")
p + theme(axis.text.x = element_text(angle=90,vjust=0.5))+ ggtitle(label=" 2024")



#find the QTLs
p <- LD.plot(datagwas2, max.loci=1000)
p + xlim(0,100) 

library(knitr)
qtl <- get.QTL(data=datagwas2,traits=c("DryWeight","SpecificGravity"),models="1-dom",bp.window=5e6)
knitr::kable(qtl)

fit.ans <- fit.QTL(data=datagwas2,trait="DryWeight",
                   qtl=qtl[,c("Marker","Model")],
                   fixed=data.frame(Effect="env",Type="factor"))
knitr::kable(fit.ans,digits=3)


