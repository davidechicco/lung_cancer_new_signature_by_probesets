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
GSE_code <- "GSE101929"
thisGEOplatform <- "GPL570"
datasetName <-  "Mitchell"
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
    label_list <- as.numeric(revalue(gset@phenoData@data$"characteristics_ch1.6", c("tumor_normal status: N"=0, "tumor_normal status: T"=1)))
    
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



    # ourSignature <- c("ADH1A",     "ANK3",      "ARF6",      "BICD2",      "FAM106A",      "FARP2",      "FKBP1A",      "HLA-DQA1",      "N4BP2L1",      "NPAS2",      "P2RX4",     "PCSK5",      "RBP1",      "RPL10L",      "SLC16A4",      "SLC39A14",      "SOX15",      "SPATA20",      "TMEM45A",      "TRIM3",      "TSFM",      "WARS2")

#     ourSignature <- c("203423_at",   "204088_at",   "204511_at",   "204911_s_at", "205234_at",  
#  "205560_at",   "206122_at",   "207820_at",   "209442_x_at", "210187_at",  
#  "212110_at",   "212656_at",   "212702_s_at", "213375_s_at", "213462_at",  
#  "213831_at",   "214182_at",   "217559_at",   "218164_at",   "218766_s_at",
#  "219410_at",   "220575_at")
 
     RotunnoSignature <- c(" 204198_s_at", "204731_at", "215806_x_at", "209813_x_at", "211144_x_at", "216920_s_at", "215227_x_at", "211571_s_at", "215646_s_at", "36936_at")
    
    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% RotunnoSignature, ])

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

source("./svm_poly_kernel.r")