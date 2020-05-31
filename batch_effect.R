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

#args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/paired-end/normalised_cpm.csv",
#          "/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary.csv",
#          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/paired-end")

args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/paired-end/processed/mode_II/all_counts_PE.csv",
          "/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv",
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

IDs_final <- sub("X*", "", colnames(df))
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
pca <- prcomp(df_ano_final, center = TRUE, scale. = FALSE)

# create final object to pass to plot
pca.final <- as.data.frame(pca$x)
pca.final$group <- ano

# make plots
# just a plot
pl <- ggplot(pca.final,aes(x=PC1,y=PC2,color=group))

# no labels
pl1 <- pl+
       geom_point(size = 3)+
       ggtitle("PCA - batch effect [paired-end]")+
       theme(plot.title = element_text(hjust = 0.5, size= 18))

# with labels
pl2 <- pl+
       geom_point(size = 3)+
       ggtitle("PCA - batch effect [paired-end]")+
       theme(plot.title = element_text(hjust = 0.5, size= 18))+
       geom_text_repel(aes(label=rownames(df_ano_final)),size=4)

# with selected labels (outliers)
pca.final$name <- ""
outliers <- c(8546, 14914, 12388, 16098, 12331, 14925, 16649)
outliers_id <- which(row.names(pca.final) %in% outliers)

pca.final$name[outliers_id] <- rownames(pca.final)[outliers_id]
pl_bis <- ggplot(pca.final,aes(x=PC1,y=PC2,color=group))
pl3 <- pl_bis+ 
       geom_point(size = 3)+
       ggtitle("PCA - batch effect [paired-end]")+
       theme(plot.title = element_text(hjust = 0.5, size= 18))+
       geom_text_repel(aes(label=name),size=4)

# write each plot as jpg file
setwd(paste(args[3], "/batch_effect", sep = ""))
jpeg('PCA_batch_effect_PE_no_labels_NEW.jpg'); pl1; dev.off()
jpeg('PCA_batch_effect_PE_labels.jpg'); pl2; dev.off()
jpeg('PCA_batch_effect_PE_outliers_labelled.jpg'); pl3; dev.off()

# remove outliers and make PCA again:
# => without 8546
# remove row with 8546 from df_ano_final and ano
to_be_removed_1 <- which(rownames(df_ano_final) == "8546")
df_ano_final_mod1 <- df_ano_final[-to_be_removed_1,]
ano_mod1 <- ano[-to_be_removed_1]

# PCA with centering and scaling
pca_mod1 <- prcomp(df_ano_final_mod1, center = TRUE, scale. = TRUE)

# create final object to pass to plot
pca.final_mod1 <- as.data.frame(pca_mod1$x)
pca.final_mod1$group <- ano_mod1
#pca.final_mod1$group <- conditionGroupsForeachBamFile

# make plots
# just a plot
pl_mod1 <- ggplot(pca.final_mod1,aes(x=PC1,y=PC2,color=group))

pl2_mod1 <- pl_mod1+
            geom_point(size = 3)+
            ggtitle("PCA - [PE] (8546 excluded)")+
            theme(plot.title = element_text(hjust = 0.5, size= 18))+
            geom_text_repel(aes(label=rownames(df_ano_final_mod1)),size=4)

# => without 8546 and 6 other: 12331, 12388, 14914, 14925, 16098, 16649

# remove row with 8546, 12331, 12388, 14914, 14925, 16098, 16649 from df_ano_final and ano
to_be_removed_2 <- which(rownames(df_ano_final) %in% c("8546", "12331", "12388", "14914", "14925", "16098", "16649"))
df_ano_final_mod2 <- df_ano_final[-to_be_removed_2,]
ano_mod2 <- ano[-to_be_removed_2]

# PCA with centering and scaling
pca_mod2 <- prcomp(df_ano_final_mod2, center = TRUE, scale. = TRUE)

# create final object to pass to plot
pca.final_mod2 <- as.data.frame(pca_mod2$x)
pca.final_mod2$group <- ano_mod2
#pca.final_mod1$group <- conditionGroupsForeachBamFile

# make plots
# just a plot
pl_mod2 <- ggplot(pca.final_mod2,aes(x=PC1,y=PC2,color=group))

pl2_mod2 <- pl_mod2+
  geom_point(size = 3)+
  ggtitle("PCA - [PE] (8546, 12331, 12388, 14914, 14925, 16098, 16649 excluded)")+
  theme(plot.title = element_text(hjust = 0.5, size= 18))+
  geom_text_repel(aes(label=rownames(df_ano_final_mod2)),size=4)

# => without 8546 and 3 other: 12388, 14914, 16098

# remove row with 8546,12388, 14914, 16098 from df_ano_final and ano
to_be_removed_3 <- which(rownames(df_ano_final) %in% c("8546", "12388", "14914", "16098"))
df_ano_final_mod3 <- df_ano_final[-to_be_removed_3,]
ano_mod3 <- ano[-to_be_removed_3]

# PCA with centering and scaling
pca_mod3 <- prcomp(df_ano_final_mod3, center = TRUE, scale. = TRUE)

# create final object to pass to plot
pca.final_mod3 <- as.data.frame(pca_mod3$x)
pca.final_mod3$group <- ano_mod3
#pca.final_mod1$group <- conditionGroupsForeachBamFile

# make plots
# just a plot
pl_mod3 <- ggplot(pca.final_mod3,aes(x=PC1,y=PC2,color=group))

pl2_mod3 <- pl_mod3+
  geom_point(size = 3)+
  ggtitle("PCA - [PE] (8546, 12388, 14914, 16098 excluded)")+
  theme(plot.title = element_text(hjust = 0.5, size= 18))+
  geom_text_repel(aes(label=rownames(df_ano_final_mod3)),size=4)

# save all the plots with outliers removed
setwd(paste(args[3], "/no_outliers", sep = ""))
jpeg('PCA_batch_effect_PE_8546_excluded.jpg'); pl2_mod1; dev.off()
jpeg('PCA_batch_effect_PE_7_samples_excluded.jpg'); pl2_mod2; dev.off()
jpeg('PCA_batch_effect_PE_4_samples_excluded.jpg'); pl2_mod3; dev.off()

################################################################################################
# label according to genders
condycje <- anno$gender..1.male..2.female.
condycje_1 <- sub(1, "male", condycje)
set_gender <- sub(2, "female", condycje_1)

# remove row with 8546,12388, 14914, 16098 from df_ano_final and ano
df_ano_final_mod4 <- df_ano_final[-to_be_removed_3,]
ano_mod4 <- set_gender[-to_be_removed_3]

# PCA with centering and scaling
pca_mod4 <- prcomp(df_ano_final_mod4, center = TRUE, scale. = TRUE)

# create final object to pass to plot
pca.final_mod4 <- as.data.frame(pca_mod4$x)
pca.final_mod4$group <- ano_mod4
#pca.final_mod1$group <- conditionGroupsForeachBamFile

# make plots
# just a plot
pl_mod4_1vs2 <- ggplot(pca.final_mod4,aes(x=PC1,y=PC2,color=group))
pl_mod4_3vs4 <- ggplot(pca.final_mod4,aes(x=PC3,y=PC4,color=group))

pl2_mod4_1vs2 <- pl_mod4_1vs2+
  geom_point(size = 3)+
  ggtitle("PCA - [PE] gender comparison")+
  theme(plot.title = element_text(hjust = 0.5, size= 18))
#  geom_text_repel(aes(label=rownames(df_ano_final_mod4)),size=4)

pl2_mod4_3vs4 <- pl_mod4_3vs4+
  geom_point(size = 3)+
  ggtitle("PCA - [PE] gender comparison")+
  theme(plot.title = element_text(hjust = 0.5, size= 18))

# variance explained plot
#By dividing the variances by the sum, we get a plot of the ratio of variance explained by each principal component.
mx_transformed <- pca_mod4$x

#These are the sample variances of the new variables.
vars_transformed <- apply(mx_transformed, 2, var)
# or: pca$sdev^2

# Notice that their sum, the total variance, is the same as for the original variables: 2.09.
#And these are the same variances divided by the total variance, i.e. how much of the total variance each new variable explains:
  
variance_explained <- vars_transformed/sum(vars_transformed)

# write plots
#variance_explained[1:4]
setwd(paste(args[3], "/gender", sep = ""))
jpeg('PCA_gender_PE_PC1vsPC2.jpg'); pl2_mod4_1vs2; dev.off()
jpeg('PCA_gender_PE_PC3vsPC4.jpg'); pl2_mod4_3vs4; dev.off()


################################################################################################
# label according to visual.loss.at.BL 
visual_loss <- anno$visual.loss.at.BL..0.no..1.yes.
visual_loss_1 <- sub(0, "no", visual_loss)
set_visual_loss <- sub(1, "yes", visual_loss_1)

# remove row with 8546,12388, 14914, 16098 from df_ano_final and ano
df_ano_final_mod5 <- df_ano_final[-to_be_removed_3,]
ano_mod5 <- set_visual_loss[-to_be_removed_3]

# PCA with centering and scaling
pca_mod5 <- prcomp(df_ano_final_mod5, center = TRUE, scale. = TRUE)

# create final object to pass to plot
pca.final_mod5 <- as.data.frame(pca_mod5$x)
pca.final_mod5$group <- ano_mod5
#pca.final_mod1$group <- conditionGroupsForeachBamFile

# make plots
# just a plot
pl_mod5_1vs2 <- ggplot(pca.final_mod5,aes(x=PC1,y=PC2,color=group))
pl_mod5_3vs4 <- ggplot(pca.final_mod5,aes(x=PC3,y=PC4,color=group))

pl2_mod5_1vs2 <- pl_mod5_1vs2+
  geom_point(size = 3)+
  ggtitle("PCA - [PE] visual loss at BL")+
  theme(plot.title = element_text(hjust = 0.5, size= 18))
#  geom_text_repel(aes(label=rownames(df_ano_final_mod4)),size=4)

pl2_mod5_3vs4 <- pl_mod5_3vs4+
  geom_point(size = 3)+
  ggtitle("PCA - [PE] visual loss at BL")+
  theme(plot.title = element_text(hjust = 0.5, size= 18))

# write plots
setwd(paste(args[3], "/visual_loss", sep = ""))
jpeg('PCA_visual_loss_PE_PC1vsPC2.jpg'); pl2_mod5_1vs2; dev.off()
jpeg('PCA_visual_loss_PE_PC3vsPC4.jpg'); pl2_mod5_3vs4; dev.off()

################################################################################################
# label according to jaw.claudication.at.BL 
jaw_claudication <- anno$jaw.claudication.at.BL...0.no..1.yes.
jaw_claudication_1 <- sub(0, "no", jaw_claudication)
set_jaw_claudication <- sub(1, "yes", jaw_claudication_1)

# remove row with 8546,12388, 14914, 16098 from df_ano_final and ano
df_ano_final_mod6 <- df_ano_final[-to_be_removed_3,]
ano_mod6 <- set_jaw_claudication[-to_be_removed_3]

# PCA with centering and scaling
pca_mod6 <- prcomp(df_ano_final_mod6, center = TRUE, scale. = TRUE)

# create final object to pass to plot
pca.final_mod6 <- as.data.frame(pca_mod6$x)
pca.final_mod6$group <- ano_mod6

# make plots
# just a plot
pl_mod6_1vs2 <- ggplot(pca.final_mod6,aes(x=PC1,y=PC2,color=group))
pl_mod6_3vs4 <- ggplot(pca.final_mod6,aes(x=PC3,y=PC4,color=group))

pl2_mod6_1vs2 <- pl_mod6_1vs2+
  geom_point(size = 3)+
  ggtitle("PCA - [PE] jaw claudication at BL")+
  theme(plot.title = element_text(hjust = 0.5, size= 18))
#  geom_text_repel(aes(label=rownames(df_ano_final_mod4)),size=4)

pl2_mod6_3vs4 <- pl_mod6_3vs4+
  geom_point(size = 3)+
  ggtitle("PCA - [PE] jaw claudication at BL")+
  theme(plot.title = element_text(hjust = 0.5, size= 18))

# write plots
setwd(paste(args[3], "/jaw_claudication", sep = ""))
jpeg('PCA_jaw_claudication_PE_PC1vsPC2.jpg'); pl2_mod6_1vs2; dev.off()
jpeg('PCA_jaw_claudication_PE_PC3vsPC4.jpg'); pl2_mod6_3vs4; dev.off()


################################################################################################
# label according to ischaemic.features.at.BL
ischaemic_features <- anno$ischaemic.features.at.BL...0.no..1.yes.
ischaemic_features_1 <- sub(0, "no", ischaemic_features)
set_ischaemic_features <- sub(1, "yes", ischaemic_features_1)

# remove row with 8546,12388, 14914, 16098 from df_ano_final and ano
df_ano_final_mod7 <- df_ano_final[-to_be_removed_3,]
ano_mod7 <- set_ischaemic_features[-to_be_removed_3]

# PCA with centering and scaling
pca_mod7 <- prcomp(df_ano_final_mod7, center = TRUE, scale. = TRUE)

# create final object to pass to plot
pca.final_mod7 <- as.data.frame(pca_mod7$x)
pca.final_mod7$group <- ano_mod6

# make plots
# just a plot
pl_mod7_1vs2 <- ggplot(pca.final_mod7,aes(x=PC1,y=PC2,color=group))
pl_mod7_3vs4 <- ggplot(pca.final_mod7,aes(x=PC3,y=PC4,color=group))

pl2_mod7_1vs2 <- pl_mod7_1vs2+
  geom_point(size = 3)+
  ggtitle("PCA - [PE] ischaemic features at BL")+
  theme(plot.title = element_text(hjust = 0.5, size= 18))
#  geom_text_repel(aes(label=rownames(df_ano_final_mod4)),size=4)

pl2_mod7_3vs4 <- pl_mod7_3vs4+
  geom_point(size = 3)+
  ggtitle("PCA - [PE] ischaemic features at BL")+
  theme(plot.title = element_text(hjust = 0.5, size= 18))

# write plots
setwd(paste(args[3], "/ischaemic_features", sep = ""))
jpeg('PCA_ischaemic_features_PE_PC1vsPC2.jpg'); pl2_mod7_1vs2; dev.off()
jpeg('PCA_ischaemic_features_PE_PC3vsPC4.jpg'); pl2_mod7_3vs4; dev.off()

