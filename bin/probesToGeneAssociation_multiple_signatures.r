
setwd(".")
options(stringsAsFactors = FALSE)

list.of.packages <- c("easypackages") # other packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("easypackages")
libraries(list.of.packages)

#BiocManager::install("annotate")
#library(annotate)

listOfBiocPackages <- c("annotate")

bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")

# check there's still something left to install
if( length(bioCpackagesNotInstalled) ) {
    biocLite(listOfBiocPackages[bioCpackagesNotInstalled])
}

library("easypackages")
libraries(list.of.packages)
libraries(listOfBiocPackages)


# retrieveSymbolsOfSignature
retrieveSymbolsOfSignature <- function(theProbesList, thePlatformAnns, keyword)
{

      #cat("probeset IDs:\n")
      thePlatformAnns <- thePlatformAnns[thePlatformAnns$"ID" %in% theProbesList,]
      # cat(platform_ann_ID_list$"Gene Symbol", sep="\n")
      #cat("\n")

      cat(keyword, ":\n")
      print(thePlatformAnns[,c("ID", "Gene Symbol")], row.names = FALSE)

      # print(cat(platform_ann_ID_list[,c("Gene Symbol")], sep="\", \""))

}

thisGEOplatform <- "GPL96"

cat("platform: ",  thisGEOplatform, "\n", sep="")
platform_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)

topSignatureProbesList <- c("201856_s_at", "201911_s_at", "202688_at", "204511_at", "204999_s_at", "207730_x_at", "208471_at", "209985_s_at", "210065_s_at", "211013_x_at", "211838_x_at", "212976_at", "213439_x_at", "213770_at", "214730_s_at", "214734_at", "215066_at", "215604_x_at", "217989_at", "218164_at", "222155_s_at", "222168_at", "37152_at")

#    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% ourSignature, ])
    
retrieveSymbolsOfSignature(topSignatureProbesList, platform_ann_df, "top signature 01")

our2ndTopSignature <- c("200785_s_at", "200984_s_at", "201481_s_at", "201798_s_at", "203013_at", "205234_at", "205560_at", "206807_s_at", "206866_at", "207820_at", "209720_s_at", "209962_at", "212702_s_at", "213228_at", "214121_x_at", "214261_s_at", "215063_x_at", "215504_x_at", "217724_at", "218164_at", "218380_at", "219410_at", "220389_at", "34206_at")
retrieveSymbolsOfSignature(our2ndTopSignature, platform_ann_df, "top signature 02")


our3rdTopSignature <- c("201481_s_at", "201729_s_at", "201798_s_at", "202418_at", "203413_at", "204999_s_at", "208082_x_at", "208471_at", "208892_s_at", "209205_s_at", "211838_x_at", "212702_s_at", "212958_x_at", "213439_x_at", "213778_x_at", "214734_at", "215066_at", "215407_s_at", "215845_x_at", "217724_at", "217989_at", "219235_s_at", "219355_at", "221668_s_at", "222168_at")
retrieveSymbolsOfSignature(our3rdTopSignature, platform_ann_df, "top signature 03")


our4thTopSignature <- c("200097_s_at", "200785_s_at", "201481_s_at", "201729_s_at", "201738_at", "202522_at", "202627_s_at", "202784_s_at", "202935_s_at", "203013_at", "203413_at", "205036_at", "207078_at", "209720_s_at", "209962_at", "213243_at", "214261_s_at", "216625_at", "217263_x_at", "217296_at", "217933_s_at", "217989_at", "218741_at", "220389_at", "222212_s_at")
retrieveSymbolsOfSignature(our4thTopSignature, platform_ann_df, "top signature 04")


our5thTopSignature <- c("200785_s_at", "201856_s_at", "201926_s_at", "202085_at", "202170_s_at", "202922_at", "204511_at", "205364_at", "206472_s_at", "207730_x_at", "208664_s_at", "208892_s_at", "210404_x_at", "211741_x_at", "212110_at", "213778_x_at", "214261_s_at", "214730_s_at", "215066_at", "216735_x_at", "218380_at", "220389_at", "221610_s_at", "222086_s_at", "222168_at")
retrieveSymbolsOfSignature(our5thTopSignature, platform_ann_df, "top signature 05")

our6thTopSignature <- c("200785_s_at", "202286_s_at", "203013_at", "203413_at", "204511_at", "205211_s_at", "206385_s_at", "206542_s_at", "210705_s_at", "210810_s_at", "211838_x_at", "213243_at", "213778_x_at", "214182_at", "215407_s_at", "217679_x_at", "217715_x_at", "217803_at", "217989_at", "218164_at", "218741_at", "219051_x_at", "219233_s_at", "220389_at", "220575_at")
retrieveSymbolsOfSignature(our6thTopSignature, platform_ann_df, "top signature 06")

