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

# BatySignatureSymbols <- c("FTSJ1", "CLIP1", "INO80", "MGC10067", "AP3D1", "MBTPS1", "VEGFB", "MYO1E", "NDUFC1", "ARG2", "NOC4", "MUS81", "MOSPD3", "HEBP2", "SDF2", "CYB561D2", "TCEB3", "MUT", "RNF103", "SNAP29", "MRPL44", "ACTR10", "OPTN", "FOXQ1", "PMS2L9", "LRRC9", "XP_372900.2", "LOC91661", "MGC33302", "SLC37A2", "MAT2B", "LOC90410", "Q6PIE2", "APG5L", "CSNK1A1", "ARPC2", "NYREN18", "CHCHD2", "RMI1", "KIAA2010", "CBWD2")

 BatySignatureSymbols <- c("ATG5", "PP4R3A", "NOC4L", "NUB1", "PMS2P3", "NPIPB3")


    # we retrieve the platform details
    platform_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
    platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
    print("sort(names(platform_ann_df))")
    print(sort(names(platform_ann_df)))

    geneSymbolFlag <- TRUE
    geneSymbolAssignment <- FALSE
    
    emptyGeneSymbol <- ""
    FIRST_GENE_EXPRESSION_INDEX <- 2
    
   thisGene <- platform_ann_df[platform_ann_df$"Gene Symbol" %in% BatySignatureSymbols,]
         
print(thisGene)

cat(thisGene[,1], sep="\n")