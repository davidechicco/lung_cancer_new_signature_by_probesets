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

ZanderSignatureData <- read.csv("../data/Zander_signature/Zander_signature_GPL6102_probesets_sequences_symbols_rand8176.csv", header = TRUE, sep =",", stringsAsFactors=FALSE)


ZanderSignatureSymbols <- ZanderSignatureData$"Symbol"




    # we retrieve the platform details
    platform_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
    platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
    print("sort(names(platform_ann_df))")
    print(sort(names(platform_ann_df)))

    geneSymbolFlag <- TRUE
    geneSymbolAssignment <- FALSE
    
    emptyGeneSymbol <- ""
    FIRST_GENE_EXPRESSION_INDEX <- 2
    
   theseGenes <- platform_ann_df[platform_ann_df$"Gene Symbol" %in% ZanderSignatureSymbols,]
         
print(theseGenes[,c("ID", "Gene Symbol")])