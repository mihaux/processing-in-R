# this script performs PCA to visulize batch effect and then remove it using removeBatchEffect() function [edgeR package]
# source: https://support.bioconductor.org/p/87718/ 

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
library(ggplot2)
library(ggrepel)
library(gridExtra)
#if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
#if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma"); library(limma)

#args <- commandArgs(trailingOnly = TRUE)

args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/paired-end/nomalised_cpm.csv",
          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/CiC_Clinical_data_FINAL.csv",
          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/paired-end")

cat("Example of usage: \n Rscript batch_effect.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end/nomalised_cpm.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/CiC_Clinical_data_FINAL.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end")

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to .csv file with normalised count data, 
       \n(2 - annotation) path to .csv annotation file (clinical data)
       \nand (3 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[1], sep="\n")

cat("Directory for results (OUT): ")
cat(args[3], sep="\n")

setwd(args[3])

# load count data
df <- read.csv(args[1], row.names = 1, header = TRUE)

IDs_final <- sub("X*", "ID_", colnames(df))
rows_final <- sub("*", "gene_", rownames(df))

colnames(df) <- IDs_final
rownames(df) <- rows_final

# the colnames need to be changed, as we are running for mode_II 
#IDs <- sub(".Aligned.sortedByCoord.out.bam*", "", colnames(df))
#IDs_final <- sub("X*", "", IDs)
#colnames(df) <- IDs_final

# load annotation (clinical) data
anno <- read.csv(args[2], row.names = 1)


# create a vector with groups based on batch information
#group_temp <- gsub("one", 1, anno$Batch)
#group <- gsub("two", 2, group_temp)

# create a vector with groups
ano <- as.vector(anno$Batch)

# add labels
df_ano <- rbind(df, ano)
rownames(df_ano)[nrow(df_ano)] <- "batch"

# omit NA values (not needed here, since there's no NAs in the dataset)
#df_ano_nona <- na.omit(df_ano)

# transpose
df_ano_tr <- data.frame(t(df_ano)) #  41 17831

# make sure your data is numeric
df_mat <- as.matrix(df_ano_tr)

# 1 - 10
df_ano_final_1 <- mapply(as.character(df_mat[1,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_2 <- mapply(as.character(df_mat[2,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_3 <- mapply(as.character(df_mat[3,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_4 <- mapply(as.character(df_mat[4,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_5 <- mapply(as.character(df_mat[5,1:(nrow(df_ano)-1)]), FUN=as.numeric)

df_ano_final_6 <- mapply(as.character(df_mat[6,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_7 <- mapply(as.character(df_mat[7,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_8 <- mapply(as.character(df_mat[8,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_9 <- mapply(as.character(df_mat[9,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_10 <- mapply(as.character(df_mat[10,1:(nrow(df_ano)-1)]), FUN=as.numeric)

# 11 - 20
df_ano_final_11 <- mapply(as.character(df_mat[11,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_12 <- mapply(as.character(df_mat[12,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_13 <- mapply(as.character(df_mat[13,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_14 <- mapply(as.character(df_mat[14,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_15 <- mapply(as.character(df_mat[15,1:(nrow(df_ano)-1)]), FUN=as.numeric)

df_ano_final_16 <- mapply(as.character(df_mat[16,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_17 <- mapply(as.character(df_mat[17,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_18 <- mapply(as.character(df_mat[18,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_19 <- mapply(as.character(df_mat[19,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_20 <- mapply(as.character(df_mat[20,1:(nrow(df_ano)-1)]), FUN=as.numeric)

# 21 - 30
df_ano_final_21 <- mapply(as.character(df_mat[21,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_22 <- mapply(as.character(df_mat[22,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_23 <- mapply(as.character(df_mat[23,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_24 <- mapply(as.character(df_mat[24,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_25 <- mapply(as.character(df_mat[25,1:(nrow(df_ano)-1)]), FUN=as.numeric)

df_ano_final_26 <- mapply(as.character(df_mat[26,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_27 <- mapply(as.character(df_mat[27,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_28 <- mapply(as.character(df_mat[28,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_29 <- mapply(as.character(df_mat[29,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_30 <- mapply(as.character(df_mat[30,1:(nrow(df_ano)-1)]), FUN=as.numeric)

# 31 - 40
df_ano_final_31 <- mapply(as.character(df_mat[31,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_32 <- mapply(as.character(df_mat[32,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_33 <- mapply(as.character(df_mat[33,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_34 <- mapply(as.character(df_mat[34,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_35 <- mapply(as.character(df_mat[35,1:(nrow(df_ano)-1)]), FUN=as.numeric)

df_ano_final_36 <- mapply(as.character(df_mat[36,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_37 <- mapply(as.character(df_mat[37,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_38 <- mapply(as.character(df_mat[38,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_39 <- mapply(as.character(df_mat[39,1:(nrow(df_ano)-1)]), FUN=as.numeric)
df_ano_final_40 <- mapply(as.character(df_mat[40,1:(nrow(df_ano)-1)]), FUN=as.numeric)

df_ano_final_41 <- mapply(as.character(df_mat[41,1:(nrow(df_ano)-1)]), FUN=as.numeric)

df_ano_final <- rbind(df_ano_final_1, df_ano_final_2, df_ano_final_3, df_ano_final_4, df_ano_final_5,
                      df_ano_final_6, df_ano_final_7, df_ano_final_8, df_ano_final_9, df_ano_final_10,
                      df_ano_final_11, df_ano_final_12, df_ano_final_13, df_ano_final_14, df_ano_final_15,
                      df_ano_final_16, df_ano_final_17, df_ano_final_18, df_ano_final_19, df_ano_final_20,
                      df_ano_final_21, df_ano_final_22, df_ano_final_23, df_ano_final_24, df_ano_final_25,
                      df_ano_final_26, df_ano_final_27, df_ano_final_28, df_ano_final_29, df_ano_final_30,
                      df_ano_final_31, df_ano_final_32, df_ano_final_33, df_ano_final_34, df_ano_final_35,
                      df_ano_final_36, df_ano_final_37, df_ano_final_38, df_ano_final_39, df_ano_final_40,
                      df_ano_final_41)

colnames(df_ano_final) <- rows_final
rownames(df_ano_final) <- IDs_final

# PCA with centering and scaling
pca <- prcomp(df_ano_final, center = TRUE, scale. = TRUE)

# create final object to pass to plot
pca.final <- as.data.frame(pca$x)
pca.final$group <- ano

# make plots
# no labels
pl1 <- ggplot(pca.final,aes(x=PC1,y=PC2,color=group))
pl1 <- pl1+
        geom_point(size = 5)+
        ggtitle("PCA - batch effect [paired-end]")+
        theme(plot.title = element_text(hjust = 0.5, size= 18))

# with labels
pl2 <- pl1+
        geom_point(size = 5)+
        ggtitle("PCA - batch effect [paired-end]")+
        theme(plot.title = element_text(hjust = 0.5, size= 18))+
        geom_text_repel(aes(label=rownames(df_ano_final)),size=3)


grid.arrange(pl1, pl2, nrow=2, ncol=1)


