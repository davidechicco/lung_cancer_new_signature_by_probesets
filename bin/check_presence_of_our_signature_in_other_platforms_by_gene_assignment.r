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

signatureName <- "our new signature"
thisSignature <- c("ZFR", "FARP1", "TNFSF10", "FARP2", "ATF5", "HPR", "ASCL1", "UPK1B", "PML", "PCDHA5", "LRRC8B", "RUNDC3A", "KSR1", "GLG1", "EXPH5", "PTPRF", "HSD17B11", "SPATA20", "SLC52A2", "AF198444", "PPARD")

theGEOplatform <- "GPL11532"
platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
print("sort(names(thisGEOplatformAnnotations))")
print(sort(names(thisGEOplatformAnnotations)))   
thisSymbolVariableName <- "gene_assignment"

theseGenes <- NULL

countDifferentGeneSymbols <- 0

for(k in 1:length(thisSignature)) {

    genesFound <- thisGEOplatformAnnotations[grepl(thisSignature[k], thisGEOplatformAnnotations[,thisSymbolVariableName], fixed=TRUE), ]$"ID"

    cat("Number of probesets found for the gene ", thisSignature[k], ": ", length(genesFound), "\n", sep="")
    
    if(length(genesFound) > 0) { countDifferentGeneSymbols <- countDifferentGeneSymbols + 1 }
    
    if(k == 1) {
	    theseGenes <- genesFound
      } else	 {
	  theseGenes <- c(theseGenes, genesFound)	
      }
}

cat("countDifferentGeneSymbols: ", countDifferentGeneSymbols, "\n", sep="")

percGeneSymbolsFound <- countDifferentGeneSymbols*100 / length(thisSignature)

cat("platform: ", theGEOplatform, "\n", sep="")

cat("signatureName: ", signatureName, "\n", sep="")
cat(countDifferentGeneSymbols, " gene symbols found out of ",  length(thisSignature), ", that is ", round(percGeneSymbolsFound,2), "%\n",  sep="")


