# this script performs PCA analysis after 'run_DESeq2.R'

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr"); library(stringr)
if (!requireNamespace("vsn", quietly = TRUE)) install.packages("vsn"); library(vsn)

library(clusterProfiler); library(readr); library(pheatmap); library(RColorBrewer); library(ggplot2)

# create a shortcut for the OneDrive directory where all files are stored
main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"      # on my mac
# main_dir <- "/Users/ummz/OneDrive - University of Leeds"                # on uni mac

#args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to .csv file with count data, 
       \n(2 - annotation) path to .csv annotation file 
       \nand (3 - output) path where output files should be stored", call.=FALSE)
}

#args <- c(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/run_1/DESeq2/Normalised_DESeq_normal_transformation_dataset_dups_run1_SE_x.Rda"),
#          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
#          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/run_1/PCA_DESeq2/"))

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3
args <- c(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/run_1/DESeq2/Normalised_DESeq_normal_transformation_dataset_nodups_run1_PE_x.Rda"),
          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/PCA_DESeq2/"))


# Example of usage: 
# Rscript run_PCA_DESeq2.R /ANALYSES/downstream/rerun_FINAL_July20/run_1/DESeq2/ /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/downstream/rerun_FINAL_July20/run_1/PCA_DESeq2/         

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

info_1 <- tail(unlist(str_split(args[1], "/")), n=1)
info_2 <- unlist(str_split(info_1, "_"))
info_3 <- paste(info_2[length(info_2)-3], info_2[length(info_2)-2], info_2[length(info_2)-1], sep = "_")

# get the file name and remove it from args[1]
dir_1 <- str_remove(args[1], info_1)

# load normalised count data
load(args[1])

# load DESeq results object
load(paste0(dir_1, "DESeq_standard_Result_dataset_", info_3, "_x.Rda"))

# load dds object
load(paste0(dir_1, "DESeqDataSet_", info_3, "_x.Rda"))

# not sure if all this code is necessary
dds <- estimateSizeFactors(dds)
normCounts <- counts(dds, normalized=TRUE)
normCounts.df <- data.frame(normCounts)
normCounts.df$REFSEQ <- rownames(normCounts.df)

sizeFactors(dds)

normCounts.df$baseMean <- res$baseMean
normCounts.df$log2FoldChange <- res$log2FoldChange
normCounts.df$lfcSE <- res$lfcSE
normCounts.df$stat <- res$stat
normCounts.df$pvalue <- res$pvalue
normCounts.df$padj <- res$padj
normCounts.df$RefSeq <- rownames(res)

if(FALSE){
try({
  names<- bitr(geneID=normCounts.df$RefSeq, fromType = "REFSEQ",toType = "SYMBOL", annoDb = "org.Hs.eg.db") #Marc1        
  #names <- bitr(geneID=normCounts.df$RefSeq, fromType = "REFSEQ", toType = "SYMBOL", OrgDb = "org.Hs.eg.db") #R Studio        
  head(names)
  merged <- merge(normCounts.df, names, by.x="RefSeq", by.y="REFSEQ", all.x=TRUE, all.y=FALSE)     
  print(head(merged))
  
  #write.csv(merged , file=paste(fileBaseName,"All_data.csv", sep="_"))
})


select <- which(res$padj < 0.01)
df <- as.data.frame(colData(dds)[, c("condition", "sample")])

normCounts <- assay(nt)[select,]

try({
  pheatmap(normCounts, cluster_rows = TRUE, show_rownames = FALSE,
           cluster_cols = TRUE, annotation_col = df, scale="row")
})
dev.copy(png, paste(fileBaseName,"DESeq_expression_heat_map_row.png", sep="_"))
dev.off()#
try({
  pheatmap(normCounts, cluster_rows = TRUE, show_rownames = FALSE,
           cluster_cols = TRUE, annotation_col = df, scale="none")
})

dev.copy(png, paste(fileBaseName,"expression_heat_map_none_Deseq2.png", sep="_"))
dev.off()#

sampleDists <- dist(t(assay(nt)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(nt$condition, nt$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy(png, paste(fileBaseName,"sample_heat_map_condition_Deseq2.png", sep="_"))
dev.off()#

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(nt$sample, nt$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy(png, paste(fileBaseName,"sample_heat_map_sample_Deseq2.png", sep="_"))
dev.off()#
}

# make PCA plots with normalisation
png(file=paste0(args[3], "PCA_plot_1_Deseq2_", info_3, "_x.png"))
plotPCA(nt, intgroup=c("condition", "sample"))
#dev.copy(png, paste0(args[3], "PCA_plot_1_Deseq2_", info_3, "_x.png"))
dev.off()

pcaData <- plotPCA(nt, intgroup=c("condition", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(file=paste0(args[3], "PCA_plot_2_Deseq2_", info_3, "_x.png"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#dev.copy(png, paste0(args[3], "PCA_plot_2_Deseq2_", info_3, "_x.png"))
dev.off()

png(file=paste0(args[3], "PCA_plot_3_Deseq2_", info_3, "_x.png"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  geom_text(aes(vjust = "inward", hjust = "inward", label = name)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#dev.copy(png, paste0(args[3], "PCA_plot_3_Deseq2_", info_3, "_x.png"))
dev.off()


if(FALSE){
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.copy(png, paste(fileBaseName,"boxplot_Deseq2.png", sep="_"))
dev.off()#

plotDispEsts(dds)
dev.copy(png, paste(fileBaseName,"Dispersion_plot_Deseq2.png", sep="_"))
dev.off()#

plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
dev.copy(png, paste(fileBaseName,"Rejected_gene_v_pvalue_cutOff_Deseq2.png", sep="_"))
dev.off()#


plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
dev.copy(png, paste(fileBaseName,"pValue_v_mean_counts_Deseq2.png", sep="_"))
dev.off()#

use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))
dev.copy(png, paste(fileBaseName,"rejected_due_to_multiple_testing_Deseq2.png", sep="_"))
dev.off()#
}

