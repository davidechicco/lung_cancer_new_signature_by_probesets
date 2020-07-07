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


# computeExecutionTime
countSignatureGeneSymbolsPresentInPlatform <- function(signatureName, thisSignature, theGEOplatformAnnotations, thisSymbolVariableName){

      theseGenes <- theGEOplatformAnnotations[theGEOplatformAnnotations[,thisSymbolVariableName] %in% thisSignature, ]
      # theseGenes <- platform_ann_df[platform_ann_df$"Gene Symbol" %in% LinSignature,]


      cat("signatureName: ", signatureName, "\n", sep="")
      print(theseGenes[,c("ID", thisSymbolVariableName)], sep="\n")

      numUniqueProbesetsInThisPlatform <- length(unique(theseGenes$"ID"))
      numUniqueGeneSymbolsInThisPlatform <- length(unique(theseGenes[, thisSymbolVariableName]))

      numGeneSymbolsSignature <- length(thisSignature)

      percGeneSymbolsFound <- numUniqueGeneSymbolsInThisPlatform *100 / numGeneSymbolsSignature

      # cat("signatureName: ", signatureName, "\n", sep="")
      cat(signatureName, " presence: ", numUniqueGeneSymbolsInThisPlatform, " unique gene symbols out of ", numGeneSymbolsSignature, " were found on this platform, that is ", round(percGeneSymbolsFound,2), "%\t", sep="")
      cat("Probesets found: ", numUniqueProbesetsInThisPlatform, "\n", sep="")

      return(theseGenes)
}

thisSignatureName <- "our new signature"
ourNewSignature <- c("ZFR", "FARP1", "TNFSF10", "FARP2", "ATF5", "HPR", "ASCL1", "UPK1B", "PML", "PCDHA5", "LRRC8B", "RUNDC3A", "KSR1", "GLG1", "EXPH5", "PTPRF", "HSD17B11", "SPATA20", "SLC52A2", "AF198444", "PPARD")

# # we retrieve the platform details
# theGEOplatform <- "GPL6097"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(platform_ann_df))")
# print(sort(names(platform_ann_df)))   
# symbolVariableName <- "Symbol"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
#
# # Result: our new signature presence: 13 unique gene symbols out of 21 were found on this platform, that is 61.9%	Probesets found: 17


# # we retrieve the platform details
# theGEOplatform <- "GPL80"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "Gene Symbol"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
# # Result: our new signature presence: 9 unique gene symbols out of 21 were found on this platform, that is 42.86%	Probesets found: 11


# # we retrieve the platform details
# theGEOplatform <- "GPL20115"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "GeneSymbol"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
# # our new signature presence: 19 unique gene symbols out of 21 were found on this platform, that is 90.48%	Probesets found: 32


# theGEOplatform <- "GPL6102"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "Symbol"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
 # our new signature presence: 18 unique gene symbols out of 21 were found on this platform, that is 85.71%	Probesets found: 27

 
#  theGEOplatform <- "GPL6884"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "Symbol"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
# # our new signature presence: 17 unique gene symbols out of 21 were found on this platform, that is 80.95%	Probesets found: 30

# theGEOplatform <- "GPL6883"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "Symbol"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
# # our new signature presence: 17 unique gene symbols out of 21 were found on this platform, that is 80.95%	Probesets found: 30

# theGEOplatform <- "GPL6480"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "GENE_SYMBOL"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
#  # our new signature presence: 19 unique gene symbols out of 21 were found on this platform, that is 90.48%	Probesets found: 38

 
# theGEOplatform <- "GPL20115"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "GeneSymbol"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
#  # our new signature presence: 19 unique gene symbols out of 21 were found on this platform, that is 90.48%	Probesets found: 32
#  
 
# theGEOplatform <- "GPL14550"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "GENE_SYMBOL"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
# # our new signature presence: 19 unique gene symbols out of 21 were found on this platform, that is 90.48%	Probesets found: 31
 
# theGEOplatform <- "GPL13497"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "GENE_SYMBOL"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
#  # our new signature presence: 19 unique gene symbols out of 21 were found on this platform, that is 90.48%	Probesets found: 31
 
# theGEOplatform <- "GPL13497"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "GENE_SYMBOL"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
#  # our new signature presence: 19 unique gene symbols out of 21 were found on this platform, that is 90.48%	Probesets found: 31

# theGEOplatform <- "GPL1293"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "Symbol"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
#  # our new signature presence: 12 unique gene symbols out of 21 were found on this platform, that is 57.14%	Probesets found: 12


# 
# theGEOplatform <- "GPL17077"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "GENE_SYMBOL"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
#  # our new signature presence: 19 unique gene symbols out of 21 were found on this platform, that is 90.48%	Probesets found: 32
 
 
#  theGEOplatform <- "GPL8300"
# platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
# thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(thisGEOplatformAnnotations))")
# print(sort(names(thisGEOplatformAnnotations)))   
# symbolVariableName <- "Gene Symbol"
# 
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
#  # our new signature presence: 15 unique gene symbols out of 21 were found on this platform, that is 71.43%	Probesets found: 26
 
 
  theGEOplatform <- "GPL6104"
platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
print("sort(names(thisGEOplatformAnnotations))")
print(sort(names(thisGEOplatformAnnotations)))   
symbolVariableName <- "Symbol"

countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ourNewSignature, thisGEOplatformAnnotations, symbolVariableName)
# our new signature presence: 18 unique gene symbols out of 21 were found on this platform, that is 85.71%	Probesets found: 24

