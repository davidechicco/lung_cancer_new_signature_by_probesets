
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

thisGEOplatform <- "GPL96"

gpl96_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
gpl96_ann_df <- as.data.frame(gpl96_ann, stringsAsFactors=FALSE)


probesList <- c("201856_s_at", "201911_s_at", "202688_at", "204511_at", "204999_s_at", "207730_x_at", "208471_at", "209985_s_at", "210065_s_at", "211013_x_at", "211838_x_at", "212976_at", "213439_x_at", "213770_at", "214730_s_at", "214734_at", "215066_at", "215604_x_at", "217989_at", "218164_at", "222155_s_at", "222168_at", "37152_at")

#    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% ourSignature, ])
    
cat("probeset IDs:\n")
gpl96_ann_ID_list <- gpl96_ann_df[gpl96_ann_df$"ID" %in% probesList,]
cat(gpl96_ann_ID_list$"Gene Symbol", sep="\n")
cat("\n")

print(gpl96_ann_ID_list[,c("ID", "Gene Symbol")], row.names = FALSE)

print(cat(gpl96_ann_ID_list[,c("Gene Symbol")], sep="\", \""))

probesWithoutSymbol <- c("207730_x_at", "215604_x_at");
annotations_without_symbols <- gpl96_ann_df[gpl96_ann_df$"ID" %in% probesWithoutSymbol,]
cat(annotations_without_symbols$"Gene Symbol", sep="\n")

