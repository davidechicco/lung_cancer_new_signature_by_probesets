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
countSignatureGeneSymbolsPresentInPlatform <- function(signatureName, thisSignature, theGEOplatformAnnotations){

      theseGenes <- theGEOplatformAnnotations[theGEOplatformAnnotations$"Gene Symbol" %in% thisSignature,]

      cat("signatureName: ", signatureName, "\n", sep="")
      print(theseGenes[,c("ID", "Gene Symbol")], sep="\n")

      numUniqueProbesetsInThisPlatform <- length(unique(theseGenes$"ID"))
      numUniqueGeneSymbolsInThisPlatform <- length(unique(theseGenes$"Gene Symbol"))

      numGeneSymbolsSignature <- length(thisSignature)

      percGeneSymbolsFound <- numUniqueGeneSymbolsInThisPlatform *100 / numGeneSymbolsSignature

      # cat("signatureName: ", signatureName, "\n", sep="")
      cat(signatureName, " presence: ", numUniqueGeneSymbolsInThisPlatform, " unique gene symbols out of ", numGeneSymbolsSignature, " were found on this platform, that is ", round(percGeneSymbolsFound,2), "%\t", sep="")
      cat("Probesets found: ", numUniqueProbesetsInThisPlatform, "\n", sep="")

      return(theseGenes)
}

theGEOplatform <- "GPL96"
# we retrieve the platform details
platform_ann <- readGEOAnn(GEOAccNum = theGEOplatform)
thisGEOplatformAnnotations <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# print("sort(names(platform_ann_df))")
# print(sort(names(platform_ann_df)))   

# thisSignatureName <- "Lin signature"
# LinSignature <- c("DDR2", "MMRN2", "C11orf80", "SLCO2B1", "RPS20P27", "OVCH1", "IRF8", "SLC26A2", "ACKR4")
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, LinSignature, thisGEOplatformAnnotations)
# 
# thisSignatureName <- "Peng signature"
# PengSignature <- c("KIF4A", "NUSAP1", "HJURP", "NEK2", "FANCI", "DTL", "UHRF1", "FEN1", "IQGAP3", "KIF20A", "TRIM59", "CENPL", "C16orf59", "UBE2C")
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, PengSignature, thisGEOplatformAnnotations)
# 
# thisSignatureName <- "Rohrbeck signature"
# RohrbeckSignature <- c("COL1A1", "KRT7", "PLS3", "TMSB10", "ARPC1B", "TUBA3", "TUBB2", "KRT15", "KRT5", "ICAM1", "ITGB2", "ITGA3", "DSG3", "DSC3", "SERPINH1", "PNMA2", "MSH3", "MYB", "RABL2B", "RABL2A")
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, RohrbeckSignature, thisGEOplatformAnnotations)

thisSignatureName <- "Sheng signature"
ShengSignature <- c("ALDH2", "ADH1B", "MAOA", "MAOB", "ABCA8",	"NR3C2", "MTHFD2", "TYMS", "FEN1", "MCM4", "MCM6", "UNG", "EME1", "BLM", "RAD51", "RAD54L", "CCNA2", "CCNB1", "ORC6", "CDC25C", "MAD2L1", "CCNE2", "CHEK1", "CDC45", "PTTG1", "CDC20", "E2F2", "CCNB2", "PLK1", "ORC1", "BUB1", "TTK", "CDK1", "ESPL1", "PKMYT1", "CDC25A", "CDC6", "CDK2", "CCNE1", "BUB1B", "E2F1", "AURKA", "SGOL1")
genesFound <- countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ShengSignature, thisGEOplatformAnnotations)

# thisSignatureName <- "Showe signature"
# ShoweSignature <- c("RSF1", "DYRK2", "YY1", "C19orf12", "THEM2", "TRIO", "MYADM", "BAIAP2", "ROGDI", "DNAJB14", "BRE", "TMEM41A", "C9orf64", "FAM110A", "PCNXL2", "REST", "C19orf62", "C13orf27", "ASCC3", "SLC1A5", "PTPLAD1", "MRE11A", "GTPBP10", "SERPINI2", "CREB1", "CCDC53", "USP48", "ZSCAN2")
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ShoweSignature, thisGEOplatformAnnotations)
# 
# thisSignatureName <- "Spira signature"
# SpiraSignature <- c("IL8", "CD55", "RGS1", "PLA2G4A", "C6", "DEFB1", "TPD52", "CD164",   "CXCL2", "SERPINA1", "FCGR3A", "TOB1", "DUSP6", "CCT2", "PPBP", "ATP6AP2", "PTK9",    "ANXA3", "FGF14", "DMD", "NELL2", "ACTR2", "CPNE3", "FOS", "SOX9", "UBXD2", "SPN",  "PPP2CA", "NUCKS1", "HDGF2", "ZC3H7B", "CCDC81", "ZNF354A", "ZNF160", "ZNF611",  "LMO4", "YWHAE", "TMED2", "DNAJC12", "RAB1A", "TRAM1", "LOC653471", "RPL35A",  "GLT28D1", "TSR1", "BACH2", "LOC153561", "COX5B", "DUOX1", "UBE2D2", "SENP6",  "FBXW12", "GTF2H3", "DCLRE1C", "SLC39A8", "SLC4A4", "TMEM47", "UBE2J1", "FXR1",   "ARL6IP5", "C1orf80", "ATP8B1", "AD7C-NTP", "STARD7", "FTO", "DKFZP434A0131", "FTL",   "FLJ14346", "PRR11", "KIAA0738", "ALMS1", "LOC152719")
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, SpiraSignature, thisGEOplatformAnnotations)
# 
# thisSignatureName <- "ZhangL signature"
# ZhangLSignature <- c("CCNI", "FGF19", "GREB1", "FRS2", "EGFR")
# countSignatureGeneSymbolsPresentInPlatform(thisSignatureName, ZhangLSignature, thisGEOplatformAnnotations)

