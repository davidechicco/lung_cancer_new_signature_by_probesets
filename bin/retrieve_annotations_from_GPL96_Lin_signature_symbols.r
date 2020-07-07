setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")

list.of.packages <- c("easypackages", "plyr") # other packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source("utils.r")

# For the printed files
num_to_return <- 1
exe_num <- sample(1:as.numeric(10000), num_to_return)


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
listOfBiocPackages <- c("annotate", "GEOquery")

bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")

# check there's still something left to install
if( length(bioCpackagesNotInstalled) ) {
    BiocManager::install(listOfBiocPackages[bioCpackagesNotInstalled])
}

library("easypackages")
libraries(list.of.packages)
libraries(listOfBiocPackages)

thisGEOplatform <- "GPL96"
# 

 LinSignature <- c("DDR2", "MMRN2", "C11orf80", "SLCO2B1", "RPS20P27", "OVCH1", "IRF8", "SLC26A2", "ACKR4")


# we retrieve the platform details
platform_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
print("sort(names(platform_ann_df))")
print(sort(names(platform_ann_df)))    
theseGenes <- platform_ann_df[platform_ann_df$"Gene Symbol" %in% LinSignature,]


print(theseGenes[,c("ID", "Gene Symbol")], sep="\n")

numUniqueProbesetsInThisPlatform <- length(unique(theseGenes$"ID"))
numUniqueGeneSymbolsInThisPlatform <- length(unique(theseGenes$"Gene Symbol"))

numGeneSymbolsSignature <- length(LinSignature)

percGeneSymbolsFound <- numUniqueGeneSymbolsInThisPlatform *100 / numGeneSymbolsSignature

cat("Presence: ", numUniqueGeneSymbolsInThisPlatform, " gene symbols out of ", numGeneSymbolsSignature, " were found on this platform for the Lin signature, that is ", round(percGeneSymbolsFound,2), "%\n", sep="")