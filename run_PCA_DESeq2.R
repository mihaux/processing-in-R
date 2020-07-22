# this script performs PCA analysis after 'run_DESeq2.R'

#Choice which normalisation method you want from the three above
#standard normalisation
load(paste(fileBaseName,"Normalised_DESeq_normal_transformation_dataset.Rda", sep="_"))
#variance stabilizing transformations
#load(paste(fileBaseName,"Normalised_DESeq_variance_stabilizing_transformations_dataset.Rda", sep="_"))
#Regularized logarithm
#load(paste(fileBaseName,"Normalised_DESeq_sRegularized_logarithm_dataset.Rda", sep="_"))

sizeFactors(dds)

#choice which result set you want
#standard results object
load(paste(fileBaseName,"DESeq_standard_Result_dataset.Rda", sep="_"))
#Adjust on false discovery rate
#load(paste(fileBaseName,"DESeq_Adjusted_for_FDR_Result_dataset.Rda", sep="_"))
#log2Folder contraction
#load(paste(fileBaseName,"DESeq_log2Folder_contraction.Rda", sep="_"))
#p Value at 0.05
#load(paste(fileBaseName,"DESeq_set_FDR_cutoff_to_05_Result_dataset.Rda", sep="_"))
#IHW DESeq_Independent_Hypothesis_Weighting_of_pvalue
#load(paste(fileBaseName,"DESeq_Independent_Hypothesis_Weighting_of_pvalue_Result_dataset.Rda", sep="_"))
#DESeq_set_analysis_using_apeglm_model
#load(paste(fileBaseName,"DESeq_set_analysis_using_apeglm_model_Result_dataset.Rda", sep="_"))

dds <- estimateSizeFactors(dds)
normCounts <- counts(dds, normalized=TRUE)
normCounts.df <- data.frame(normCounts)
#head(normCounts.df)
normCounts.df$REFSEQ <- rownames(normCounts.df)
#head(normCounts.df)

sizeFactors(dds)

#head(res)
normCounts.df$baseMean <- res$baseMean
normCounts.df$log2FoldChange <- res$log2FoldChange
normCounts.df$lfcSE <- res$lfcSE
normCounts.df$stat <- res$stat
normCounts.df$pvalue <- res$pvalue
normCounts.df$padj <- res$padj
normCounts.df$RefSeq <- rownames(res)
#head(normCounts.df)
try({
  names<- bitr(geneID=normCounts.df$RefSeq, fromType = "REFSEQ",toType = "SYMBOL", annoDb = "org.Hs.eg.db") #Marc1        
  #names <- bitr(geneID=normCounts.df$RefSeq, fromType = "REFSEQ", toType = "SYMBOL", OrgDb = "org.Hs.eg.db") #R Studio        
  head(names)
  merged <- merge(normCounts.df, names, by.x="RefSeq", by.y="REFSEQ", all.x=TRUE, all.y=FALSE)     
  print(head(merged))
  
  write.csv(merged , file=paste(fileBaseName,"All_data.csv", sep="_"))
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

plotPCA(nt, intgroup=c("condition", "sample"))
dev.copy(png, paste(fileBaseName,"PCA_plot_1_Deseq2.png", sep="_"))
dev.off()#

pcaData <- plotPCA(nt, intgroup=c("condition", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.copy(png, paste(fileBaseName,"PCA_plot_2_Deseq2.png", sep="_"))
dev.off()#

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  geom_text(aes(vjust = "inward", hjust = "inward", label = name)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.copy(png, paste(fileBaseName,"PCA_plot_2_Deseq2.png", sep="_"))
dev.off()#

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


#W <- res$stat
#maxCooks <- apply(assays(dds)[["cooks"]],1,max)
#idx <- !is.na(W)
#plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
#     ylab="maximum Cook's distance per gene",
#     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
#m <- ncol(dds)
#p <- 3
#abline(h=qf(.99, p, m - p))
#dev.copy(png, paste(fileBaseName,"outlier_detection_Deseq2.png", sep="_"))
#dev.off()#

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


