# script to compare results from Mann-Whitney for outliers vs. non-outliers and PCA loadings

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)
library(ggplot2)

# get working directory to recognise the machine
w_dir <- getwd()

# create a shortcut for the OneDrive directory where all files are stored
if(startsWith(w_dir, "/Users/michal")){           
  main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"    # on my mac
} else if (startsWith(w_dir, "/Users/ummz")) {    
  main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds"      # on uni mac    
} else {
  print("Unrecognised machine.")
}

setwd(paste0(main_dir, "/ANALYSES_archived/oct_nov/nov20_outliers_testing/comparison_with_PC_loadings"))

# load tables with results (p-values, p-adjusted)
DE_res_raw <- read.csv(paste0(main_dir, "/ANALYSES_archived/oct_nov/nov20_outliers_testing/raw/table_sorted_padjusted_raw_mod.csv"), row.names = 1)
DE_res_vst <- read.csv(paste0(main_dir, "/ANALYSES_archived/oct_nov/nov20_outliers_testing/vst/table_sorted_padjusted_vst_mod.csv"), row.names = 1)
DE_res_rlog <- read.csv(paste0(main_dir, "/ANALYSES_archived/oct_nov/nov20_outliers_testing/rlog/table_sorted_padjusted_rlog_mod.csv"), row.names = 1)

# get gene names of the results with p-adj < 0.05
DE_genes_raw <- rownames(DE_res_raw)[which(DE_res_raw$fdr.pvalue < 0.05)]
DE_genes_vst <- rownames(DE_res_vst)[which(DE_res_vst$fdr.pvalue < 0.05)]
DE_genes_rlog <- rownames(DE_res_rlog)[which(DE_res_rlog$fdr.pvalue < 0.05)]

# run PCA using the built-in method and data with gene symbols
# define directory with data (INPUT)
data_dir_loadings <- paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/all_chr/INPUT_counts")

# load data RAW | VST | rlog (run for one data type at a time)
load(paste0(data_dir_loadings, "/Raw_DESeq_dataset_all.Rda")); dds <- dds_all
#load(paste0(data_dir_loadings,"/Normalised_DESeq_vst_dataset_all.Rda")); dds <- vst_all
#load(paste0(data_dir_loadings,"/Normalised_DESeq_rlog_dataset_all.Rda")); dds <- rlog_all

# run PCA using the built-in function
pca_dds <- prcomp(t(assay(dds)))

# retrieve PCA loadings
loadings_dds <- as.data.frame(pca_dds$rotation)

# get PC1 and PC2 values for DE_genes_raw | DE_genes_vst | DE_genes_rlog
loadings_dds$PC1[which(rownames(loadings_dds) %in% DE_genes_raw)]
loadings_dds$PC2[which(rownames(loadings_dds) %in% DE_genes_raw)]

# make histograms of PC1 and PC2
hist(loadings_dds$PC1, breaks = 5, main = "Histogram of PC1 loadings - raw")
hist(loadings_dds$PC2, breaks = 5, main = "Histogram of PC2 loadings - raw")

# min and max value of PC1
#min(loadings_dds$PC1); max(loadings_dds$PC1)

# plot loadings of PC1
plot(loadings_dds$PC1, ylab="PC1 loadings")

# check the "extreme values"
plot(loadings_dds$PC1[which(loadings_dds$PC1>0.05)], main="Extreme values of PC1 loadings", ylab="PC1 loadings > 0.05")
text(loadings_dds$PC1[which(loadings_dds$PC1>0.05)], labels=rownames(loadings_dds)[which(loadings_dds$PC1>0.05)], cex=0.6, pos=3, col="red")   # print labels

# plot loadings of PC2
plot(loadings_dds$PC2, ylab="PC2 loadings")

plot(loadings_dds$PC2[which(loadings_dds$PC2>0.015)], main="Extreme values of PC2 loadings (positive)", ylab="PC2 loadings > 0.015")
text(loadings_dds$PC2[which(loadings_dds$PC2>0.015)], labels=rownames(loadings_dds)[which(loadings_dds$PC2>0.015)], cex=0.6, pos=3, col="red")   # print labels

plot(loadings_dds$PC2[which(loadings_dds$PC2<=-0.05 & loadings_dds$PC2>-0.8)], main="Extreme values of PC2 loadings (negative)", ylab="PC2 loadings < -0.05 & > -0.08")
text(loadings_dds$PC2[which(loadings_dds$PC2<=-0.05 & loadings_dds$PC2>-0.8)], labels=rownames(loadings_dds)[which(loadings_dds$PC2<=-0.05 & loadings_dds$PC2>-0.8)], cex=0.6, pos=3, col="red")   # print labels

# get loadings of DE genes from Mann-Whitney
which(rownames(loadings_dds) %in% DE_genes_raw)

rownames(loadings_dds)[which(rownames(loadings_dds) %in% DE_genes_raw)]

plot(loadings_dds$PC1[which(rownames(loadings_dds) %in% DE_genes_raw)], main="PC1 loadings of 30 DE genes - raw data", ylab="PC1 loadings")
plot(loadings_dds$PC2[which(rownames(loadings_dds) %in% DE_genes_raw)], main="PC2 loadings of 30 DE genes - raw data", ylab="PC2 loadings")

df_raw_final <- data.frame(genes=rownames(loadings_dds)[which(rownames(loadings_dds) %in% DE_genes_raw)],
                           PC1=loadings_dds$PC1[which(rownames(loadings_dds) %in% DE_genes_raw)],
                           PC2=loadings_dds$PC2[which(rownames(loadings_dds) %in% DE_genes_raw)])

write.csv(df_raw_final, file="df_raw_final.csv")
