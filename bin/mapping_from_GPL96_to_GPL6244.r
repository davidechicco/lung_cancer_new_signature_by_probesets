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

# let's read the conversion file
mapping_table <- read.csv(file = "../data/probeset_mapping_files/U133PlusVsHuGene_BestMatch.csv", header=T)

ourPlatformIDs <- mapping_table[,c("Name")]
otherPlatformIDs <- mapping_table[,c("Cluster")]

cat(head(ourPlatformIDs), "\n")
cat(head(otherPlatformIDs), "\n")

ourNewSignatureProbesets <-  c("201856_s_at", "201911_s_at", "202688_at", "204511_at", "204999_s_at", "207730_x_at", "208471_at", "209985_s_at", "210065_s_at", "211013_x_at", "211838_x_at", "212976_at", "213439_x_at", "213770_at", "214730_s_at", "214734_at", "215066_at", "215604_x_at", "217989_at", "218164_at", "222155_s_at", "222168_at", "37152_at")

numOurProbesetsPresent <- sum(ourNewSignatureProbesets %in% ourPlatformIDs)

cat("There are ", numOurProbesetsPresent, " probesets present out of ", length(ourNewSignatureProbesets), " which means ", round(numOurProbesetsPresent * 100 / length(ourNewSignatureProbesets),2), "%\n", sep="")
