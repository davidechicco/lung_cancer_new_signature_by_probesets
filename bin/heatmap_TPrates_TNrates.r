setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
# set.seed(11)
options(repos = list(CRAN="http://cran.rstudio.com/"))

# update R packages
list.of.packages <- c("easypackages", "ggplot2", "viridis", "hrbrthemes", "ztable", "devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library("easypackages")
libraries(list.of.packages)

devtools::install_github("cardiomoon/ztable")


# read the file
fileName <- "../results_our_other_top_signatures/recap_results_6_signatures_TPrate_TNrate.csv"
table_TPrates_TNrates <- read.table(file=fileName,head=TRUE,sep=",", row.names = 1, stringsAsFactors=FALSE)

table_TPrates <- table_TPrates_TNrates[1:6,]
table_TNrates <- table_TPrates_TNrates[7:12,]

# generate the heatmap
dataframe_TPrates <- as.data.frame(table_TPrates, stringsAsFactors = FALSE)
dataframe_TNrates <- as.data.frame(table_TNrates, stringsAsFactors = FALSE)

stackedDF_TPrates <- stack(table_TPrates)
stackedDF_TPrates$"sign" <- rownames(table_TPrates) # c("sign01", "sign02", "sign03", "sign04", "sign05", "sign06")
colnames(stackedDF_TPrates)[1]<- "TPrate"

stackedDF_TNrates <- stack(table_TNrates)
stackedDF_TNrates$"sign" <- rownames(table_TNrates) # c("sign01", "sign02", "sign03", "sign04", "sign05", "sign06")
colnames(stackedDF_TNrates)[1] <- "TNrate"

# TP rate heatmap
plot_title_TP_rates <- "TP rate (sensitivity, recall) obtained on the patients with lung cancer"
heatmap_plot_TP_rates <- ggplot(stackedDF_TPrates, aes(x=sign, y=ind, fill=TPrate)) + geom_tile() 
heatmap_plot_TP_rates <- heatmap_plot_TP_rates + scale_fill_gradientn(colours = c("white", "red", "blue"), breaks=c(0, 0.1, 0.25, 0.5, 0.75, 1), values = c(0,0.1,1)) 
heatmap_plot_TP_rates <- heatmap_plot_TP_rates + scale_y_discrete(limits = rev(levels(as.factor(stackedDF_TPrates$ind))))
heatmap_plot_TP_rates <- heatmap_plot_TP_rates + labs(y="datasets", x = "signatures")
heatmap_plot_TP_rates <- heatmap_plot_TP_rates + ggtitle(plot_title_TP_rates)
print(heatmap_plot_TP_rates)

heatmap_width <- 25
heatmap_height <- 8

heatmap_TP_rates <- paste0("../results_our_other_top_signatures/heatmap_our6signatures_TPrates_rand", exe_num)

SAVE_HEATMAP_TO_PNG_FILE_TPR <- TRUE

if(SAVE_HEATMAP_TO_PNG_FILE_TPR) {
    heatmap_file_png <- paste0(heatmap_TP_rates, ".png")
    ggsave(plot=heatmap_plot_TP_rates, filename=heatmap_file_png, width=heatmap_width, height=heatmap_height, units="cm")
    cat("Saved file ", heatmap_file_png, "\n")
  }

# TN rate heatmap
plot_title_TN_rates <- "TN rate (specificity) obtained on the healthy controls"
heatmap_plot_TN_rates <- ggplot(stackedDF_TNrates, aes(x=sign, y=ind, fill=TNrate)) + geom_tile() 
heatmap_plot_TN_rates <- heatmap_plot_TN_rates + scale_fill_gradientn(colours = c("white", "red", "blue"), breaks=c(0, 0.1, 0.25, 0.5, 0.75, 1), values = c(0,0.1,1)) 
heatmap_plot_TN_rates <- heatmap_plot_TN_rates + scale_y_discrete(limits = rev(levels(as.factor(stackedDF_TNrates$ind))))
heatmap_plot_TN_rates <- heatmap_plot_TN_rates + labs(y="datasets", x = "signatures")
heatmap_plot_TN_rates <- heatmap_plot_TN_rates + ggtitle(plot_title_TN_rates)
print(heatmap_plot_TN_rates)


heatmap_TN_rates <- paste0("../results_our_other_top_signatures/heatmap_our6signatures_TNrates_rand", exe_num)

SAVE_HEATMAP_TO_PNG_FILE_TNR <- TRUE

if(SAVE_HEATMAP_TO_PNG_FILE_TNR) {
    heatmap_file_png <- paste0(heatmap_TN_rates, ".png")
    ggsave(plot=heatmap_plot_TN_rates, filename=heatmap_file_png, width=heatmap_width, height=heatmap_height, units="cm")
    cat("Saved file ", heatmap_file_png, "\n")
  }

  
# 
# # first method to create the heatmap
# ztable(dataframe_TPrates) %>% 
#     makeHeatmap(palette="Blues") %>%
#     print(caption="Heatmap")
#     
#     
#     # first method to create the heatmap
# ztable(dataframe_TNrates) %>% 
#     makeHeatmap(palette="Blues") %>%
#     print(caption="Heatmap")