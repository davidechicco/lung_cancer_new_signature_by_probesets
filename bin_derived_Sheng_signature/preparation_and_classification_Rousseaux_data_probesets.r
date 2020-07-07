setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
# set.seed(11)

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

targetName <- "lung_cancer"
GSE_code <- "GSE30219"
thisGEOplatform <- "GPL570"
datasetName <-  "Rousseaux"
# folderPath <- paste0("../data/datasets_", thisGEOplatform, "/")
# createFolderCommand <- paste0("mkdir -p ", folderPath, "; \n", sep="")
# cat(createFolderCommand)
# system(createFolderCommand)


gset <- getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)

if (length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gene_expression <- as.data.frame(exprs(gset))

# cat("str(gset@phenoData@data)\n")
# print(str(gset@phenoData@data))
# cat("str(gset@phenoData@data)\n")

LABEL_DETECTED <- TRUE
 
 if(LABEL_DETECTED == TRUE) {
 
    # # # # we add the labels
    library("plyr")
    label_list <- as.numeric(revalue(gset@phenoData@data$source_name_ch1, c("Lung Tumour"=0, "Non Tumoral Lung"=1)))
    
#     label_list <- c()
#     i <- 1
#     for(thisTitle in gset@phenoData@data$"source_name_ch1") {
#       
# 	if(grepl("normal", thisTitle)) {
# 	      label_list[[i]] <- 0
# 	  } else {
# 	    label_list[[i]] <- 1
# 	  }
# 	    i <- i + 1
#       }
    
    cat("label_list:\n")
    print(label_list)
    
    labels_df_temp <- as.data.frame(label_list)
    
    labels_df <- as.data.frame(t(labels_df_temp))
    colnames(labels_df) <- colnames(gene_expression)
    rownames(labels_df) <- targetName
    gene_expression_with_labels <- rbind(labels_df, gene_expression)
    
    gene_expression_with_labels$ID <- rownames(gene_expression_with_labels)

#     # we retrieve the platform details
#     platform_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
#     platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
#     print("sort(names(platform_ann_df))")
#     print(sort(names(platform_ann_df)))

     derivedShengSignature <- c("201425_at", "201761_at", "201930_at", "202240_at", "202330_s_at", "202589_at", "202705_at", "202870_s_at", "2028_s_at", "203213_at", "203214_x_at", "203362_s_at", "203418_at", "203554_x_at", "203755_at", "203967_at", "203968_s_at", "204041_at", "204092_s_at", "204126_s_at", "204252_at", "204267_x_at", "204388_s_at", "204389_at", "204558_at", "204695_at", "204696_s_at", "204719_at", "204767_s_at", "204768_s_at", "204817_at", "204822_at", "204947_at", "205023_at", "205024_s_at", "205034_at", "205085_at", "205167_s_at", "205259_at", "205393_s_at", "205394_at", "205733_at", "207042_at", "208079_s_at", "208080_at", "209612_s_at", "209613_s_at", "209614_at", "209642_at", "210559_s_at", "211803_at", "211804_s_at", "211814_s_at", "212141_at", "212142_at", "212741_at", "213226_at", "213523_at", "214710_s_at", "215508_at", "215509_s_at", "216275_at", "216277_at", "216914_at", "217010_s_at", "217684_at", "219105_x_at", "222036_s_at", "222037_at", "38158_at")
    
    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% derivedShengSignature, ])
   

    print("dim(patients_data_filtered_ourSignature)\n")
    print(dim(patients_data_filtered_ourSignature))

    patients_data_filtered_ourSignature <- rbind(gene_expression_with_labels[1,], patients_data_filtered_ourSignature)

    # outputFileName <- paste0(folderPath,"patients_data_filtered_ourSignature_", exe_num, ".csv")
    # write.table(patients_data_filtered_ourSignature, file=outputFileName, row.names=TRUE, col.names=TRUE, sep=",")
    # cat("saved file  ", outputFileName, "\n")

    patients_data_filtered_ourSignature_t <- as.data.frame(t(patients_data_filtered_ourSignature))

    library("dplyr")
    patients_data_filtered_ourSignature_t <- patients_data_filtered_ourSignature_t %>% select(-targetName,targetName)

    patients_data_filtered_ourSignature_t <-  slice(patients_data_filtered_ourSignature_t, 1:(n()-1))
    patients_data_filtered_ourSignature_t[,targetName] <- as.numeric(patients_data_filtered_ourSignature_t[,targetName]) 

    tableOnlyOnesAndZeros <-  rbind(patients_data_filtered_ourSignature_t[patients_data_filtered_ourSignature_t$lung_cancer==0,], patients_data_filtered_ourSignature_t[patients_data_filtered_ourSignature_t$lung_cancer==1,])

#     outputFileName <- paste0(folderPath, datasetName, "_", GSE_code, "_", thisGEOplatform, "_patients_data_filtered_ourSignature_transpose_", exe_num, ".csv")
#     write.table(tableOnlyOnesAndZeros, file=outputFileName, row.names=FALSE, col.names=TRUE, sep=",")
#     cat("saved file  ", outputFileName, "\n")
    

    patients_data <- tableOnlyOnesAndZeros

}

ROSE_OVERSAMPLE <- FALSE
source("./svm_poly_kernel.r")
