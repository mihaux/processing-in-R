# this script performs DESeq2 analysis on read counts data (output from featureCounts)
# created based DESeq2 vignette on https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# DESEq2 manual:        https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf 
# guide for beginners:  https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf

# excellent source about running DE analysis
# https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html

# nice presentation
# http://people.duke.edu/~ccc14/duke-hts-2018/_downloads/stat-GLM-model-RNA-Seq-handout.pdf

# Why do we use the negative binomial distribution for analysing RNAseq data?
# http://bridgeslab.sph.umich.edu/posts/why-do-we-use-the-negative-binomial-distribution-for-rnaseq

# The package DESeq2 provides methods to test for differential expression by use of negative binomial generalized linear models; 
# the estimates of dispersion and logarithmic fold changes incorporate data-driven prior distributions.

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr"); library(stringr)
if (!requireNamespace("vsn", quietly = TRUE)) BiocManager::install("vsn"); library(vsn)

# if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
# if (!requireNamespace("vsn", quietly = TRUE)) install.packages("vsn"); library(vsn)

# get working directory to recognise the machine
w_dir <- getwd()

# create a shortcut for the OneDrive directory where all files are stored
if(startsWith(w_dir, "/Users/michal")){           
  main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"    # on my mac
} else if (startsWith(w_dir, "/Users/ummz")) {    
  main_dir <- "/Users/ummz/OneDrive - University of Leeds"                # on uni mac    
} else {
  print("Unrecognised machine.")
}

#args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to _x.csv file with count data, 
       \n(2 - annotation) path to _x.csv annotation file 
       \nand (3 - output) path where output files should be stored", call.=FALSE)
}

#args <- c(paste0(main_dir, "/ANALYSES/rerun_FINAL_July20/run_1/featCounts_PE/all_counts_dups_run1_PE_mod.csv"),
#          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
#          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/run_1/FINAL/"))

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3

args <- c(paste0(main_dir, "/ANALYSES/comparison_with_Ian_results/rerun_5/featCounts/all_counts_dups_rr5_x.csv"),
          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL/"))

# Example of usage: 
# Rscript run_DESeq2.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_x.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED_x.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/downstream/rerun_FINAL_July20/run_1/DESeq2/            

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load count data (cols -> samples | rows -> genes)
df <- read.csv(args[1], row.names = 1, header = TRUE)     # data.frame with counts only

# NOTICE: It is important to never provide counts that were pre-normalized for sequencing depth/library size, 
# as the statistical model is most powerful when applied to un-normalized counts, 
# and is designed to account for library size differences internally.

# NOTICE: the counts table must be in the form of a matrix of integer values
# switch from data.frame to matrix and to integer
counts_mat <- as.matrix(df)                       # class: matrix | type: double  | dim: 28278    41
storage.mode(counts_mat) <- "integer"             # class: matrix | type: integer

# load annotation (clinical) data which will be used for coldata
anno <- read.csv(args[2], row.names = 1)

# add "ID_" to all rownames
rownames(anno) <- paste0("ID_", rownames(anno))

# make df colnames matchin the annotation rownames
colnames(counts_mat) <- str_replace(colnames(df), "X", "ID_")

# create DESeqDataSet object from Count matrix input
# NOTE: A DESeqDataSet object must have an associated design formula. 
# The design formula expresses the variables which will be used in modeling (and is used to estimate the dispersions and to estimate the log2 fold changes of the model.)

# NOTE: It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. 

# make sure the order is the same for counts_mat and anno and names correspond
all(rownames(anno) %in% colnames(counts_mat))
all(rownames(anno) == colnames(counts_mat))

coldata <- cbind(rownames(anno), anno[,c("visual.loss.ever...0.no..1.yes." ,"gender..1.male..2.female.")])

names(coldata) <- c("sample", "visual_loss", "gender")

# transform coldata feature to text (not sure if it's needed)
coldata$visual_loss <- str_replace(coldata$visual_loss, pattern = "0", "no")
coldata$visual_loss <- str_replace(coldata$visual_loss, pattern = "1", "yes")

coldata$gender <- str_replace(coldata$gender, pattern = "1", "male")
coldata$gender <-str_replace(coldata$gender, pattern = "2", "female")

coldata$visual_loss <- factor(coldata$visual_loss)
coldata$gender <- factor(coldata$gender)


#################################################
# split df in 
# => chr 1-22
# => chr X and Y
# => all chr (no chenges needed)

# load lists with positions on chrX and chrY
chrX_IDs <- read.csv(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/PCA_plots/chrX_positions.csv"), row.names = 1)
chrY_IDs <- read.csv(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/PCA_plots/chrY_positions.csv"), row.names = 1)

# get union of chrX (l=2276) and chrY (l=303)
both_chr_IDs <- union(chrX_IDs$x, chrY_IDs$x) # l=2516

# subset df for chrX and chrY only
df_chrXY <- df[which(rownames(df) %in% both_chr_IDs),]
counts_mat_chrXY <- as.matrix(df_chrXY)                       
storage.mode(counts_mat_chrXY) <- "integer"                   # type: integer | dim: 2516    41

# subset df for chr1-22 only
df_chr1_22 <- df[-which(rownames(df) %in% both_chr_IDs),]
counts_mat_chr1_22 <- as.matrix(df_chr1_22)                   
storage.mode(counts_mat_chr1_22) <- "integer"                 # type: integer | dim: 56092    41

# df includes all chromosomes

#################################################

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~ gender)
dds_chrXY <- DESeqDataSetFromMatrix(countData = counts_mat_chrXY, colData = coldata, design = ~ gender)
dds_chr1_22 <- DESeqDataSetFromMatrix(countData = counts_mat_chr1_22, colData = coldata, design = ~ gender)

# adding additional features (if needed) to the metadata columns
#featureData <- data.frame(gene=rownames(cts))
#mcols(dds) <- DataFrame(mcols(dds), featureData)
#mcols(dds)

### pre-filtering ###
# NOTE: it's not necessary to pre-filter low count genes before running the DESeq2 functions, 
# but by removing rows in which there are very few reads, 
# we reduce the memory size of the dds data object, and 
# we increase the speed of the transformation and testing functions within DESeq2. 

# Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total. 
# Note that more strict filtering to increase power is automatically applied via independent filtering on the mean of normalized counts within the results function.

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]                                       # 58608 - 52239 = 6369 rows dropped

keep_chrXY <- rowSums(counts(dds_chrXY)) >= 10
dds_chrXY <- dds_chrXY[keep_chrXY,]                     # 2516 - 2249 = 267 rows dropped

keep_chr1_22 <- rowSums(counts(dds_chr1_22)) >= 10
dds_chr1_22 <- dds_chr1_22[keep_chr1_22,]               # 56092 - 49990 = 6102 rows dropped

# set "female" as reference level (no reason for this)
dds$gender <- factor(dds$gender, levels = c("female","male"))
dds_chrXY$gender <- factor(ddss_chrXY$gender, levels = c("female","male"))
dds_chr1_22$gender <- factor(dds_chr1_22$gender, levels = c("female","male"))

# PCAplot to see if they cluster togther
# TODO

# perform DESeq analysis
dds <- DESeq(dds)
res <- results(dds)       # dim = 52239     6

dds_chrXY <- DESeq(dds_chrXY)
res_chrXY <- results(dds_chrXY)       # dim = 52239     6

dds_chr1_22 <- DESeq(dds_chr1_22)
res_chr1_22 <- results(dds_chr1_22)       # dim = 52239     6

# resultsNames(dds)  # [1] "Intercept"             "gender_male_vs_female"

#mcols(res, use.names=TRUE)

#DataFrame with 6 rows and 2 columns
# type                                          description
#                     <character>                                   <character>
#   baseMean          intermediate        mean of normalized counts for all samples
#   log2FoldChange         results    log2 fold change (MLE): gender male vs female
#   lfcSE                  results            standard error: gender male vs female
#   stat                   results            Wald statistic: gender male vs female
#   pvalue                 results         Wald test p-value: gender male vs female
#   padj                   results                             BH adjusted p-values


# retrieve unnormalised counts they are the same before and after running DESeq()
counts_non_norm <- counts(dds)
write.csv(counts_non_norm, file="non_norm_counts.csv")

# retrive normalized counts
counts_norm <- counts(dds, normalized = TRUE)
write.csv(counts_norm, file="norm_counts.csv")

counts_norm_chrXY <- counts(dds_chrXY, normalized = TRUE)
write.csv(counts_norm_chrXY, file="norm_counts_chrXY.csv")

counts_norm_chr1_22 <- counts(dds_chr1_22, normalized = TRUE)
write.csv(counts_norm_chr1_22, file="norm_counts_chr1_22.csv")

# NOTE: PCA function needs a DESeqTransform class object => DESeqTransform(SummarizedExperiment)

# unnormalised
dds_all <- dds
dds_chr1_22
dds_chrXY

# it takes too long, especially for rlog, so just load saved data (below)
# vst (variance stabilizing transformation)
#vst_all     <- vst(dds, blind=FALSE)
#vst_chr1_22 <- vst(dds_chr1_22, blind=FALSE)
#vst_chrXY   <- vst(dds_chrXY, blind=FALSE)

# rlog (it takes quite long time to run)
#rlog_all       <- rlog(dds, blind=FALSE)
#rlog_chr1_22   <- rlog(dds_chr1_22, blind=FALSE)
#rlog_chrXY     <- rlog(dds_chrXY, blind=FALSE)

# save transformed data and unnormalised as well
save(dds_all, file="Raw_DESeq_dataset_all.Rda")
save(dds_chr1_22, file="Raw_DESeq_dataset_chr1_22.Rda")
save(dds_chrXY, file="Raw_DESeq_dataset_chrXY.Rda")

#save(vst_all, file="Normalised_DESeq_vst_dataset_all.Rda")
#save(vst_chr1_22, file="Normalised_DESeq_vst_dataset_chr1_22.Rda")
#save(vst_chrXY, file="Normalised_DESeq_vst_dataset_chrXY.Rda")

#save(rlog_all, file="Normalised_DESeq_rlog_dataset_all.Rda")
#save(rlog_chr1_22, file= "Normalised_DESeq_rlog_dataset_chr1_22.Rda")
#save(rlog_chrXY, file="Normalised_DESeq_rlog_dataset_chrXY.Rda")

# load transformed data
load(paste0(main_dir,"/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL/Normalised_DESeq_vst_dataset_all.Rda"))
load(paste0(main_dir,"/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL/Normalised_DESeq_vst_dataset_chr1_22.Rda"))
load(paste0(main_dir,"/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL/Normalised_DESeq_vst_dataset_chrXY.Rda"))

load(paste0(main_dir,"/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL/Normalised_DESeq_rlog_dataset_all.Rda"))
load(paste0(main_dir,"/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL/Normalised_DESeq_rlog_dataset_chr1_22.Rda"))
load(paste0(main_dir,"/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL/Normalised_DESeq_rlog_dataset_chrXY.Rda"))

# all the objects need to be passed through DESeqTransform() function before PCAplot()
dds_all_trans       <- DESeqTransform(dds_all)
dds_chrXY_trans     <- DESeqTransform(dds_chrXY)
dds_chr1_22_trans   <- DESeqTransform(dds_chr1_22)

vst_all_trans       <- DESeqTransform(vst_all)
vst_chr1_22_trans   <- DESeqTransform(vst_chr1_22)
vst_chrXY_trans     <- DESeqTransform(vst_chrXY)

rlog_all_trans       <- DESeqTransform(rlog_all)
rlog_chr1_22_trans   <- DESeqTransform(rlog_chr1_22)
rlog_chrXY_trans     <- DESeqTransform(rlog_chrXY)

# PCA plots
p1 <- plotPCA(dds_all_trans, intgroup=c("gender")) + ggtitle("raw: all chromosomes") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), aspect.ratio=1)
p2 <- plotPCA(dds_chrXY_trans, intgroup=c("gender")) + ggtitle("raw: chrX & chrY") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), aspect.ratio=1)
p3 <- plotPCA(dds_chr1_22_trans, intgroup=c("gender")) + ggtitle("raw: chr1 - chr22") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), aspect.ratio=1)

# vst => variance stabilizing transformation
p4 <- plotPCA(vst_all_trans, intgroup=c("gender")) + ggtitle("vst: all chromosomes") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), aspect.ratio=1)
p5 <- plotPCA(vst_chrXY_trans, intgroup=c("gender")) + ggtitle("vst: chrX & chrY") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), aspect.ratio=1)
p6 <- plotPCA(vst_chr1_22_trans, intgroup=c("gender")) + ggtitle("vst: chr1 - chr22") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), aspect.ratio=1)
  
p7 <- plotPCA(rlog_all_trans, intgroup=c("gender")) + ggtitle("rlog: all chromosomes") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1)
p8 <- plotPCA(rlog_chrXY_trans, intgroup=c("gender")) + ggtitle("rlog: chrX - chrY") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1)
p9 <- plotPCA(rlog_chr1_22_trans, intgroup=c("gender")) + ggtitle("rlog: chr1 - chr22") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1)

# plot grouped:

# unnormalised only
grid.arrange(p1, p2, p3, nrow = 1)

# normalised (vst + rlog)
grid.arrange(p4, p5, p6, p7, p8, p9, nrow = 2)

# all
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)

# row and vst only
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

# plots without titles
p1_bis <- plotPCA(dds_all_trans, intgroup=c("gender")) + theme(legend.position = "none")
p2_bis <- plotPCA(dds_chrXY_trans, intgroup=c("gender")) + theme(legend.position = "none")
p3_bis <- plotPCA(dds_chr1_22_trans, intgroup=c("gender")) + theme(legend.position = "none")

# vst => variance stabilizing transformation
p4_bis <- plotPCA(vst_all_trans, intgroup=c("gender")) + theme(legend.position = "none")
p5_bis <- plotPCA(vst_chrXY_trans, intgroup=c("gender")) + theme(legend.position = "none")
p6_bis <- plotPCA(vst_chr1_22_trans, intgroup=c("gender")) + theme(legend.position = "none")

# save plots
pdf(file="PCA_plot_gender_raw_all_chr.pdf") ; p1; dev.off()
pdf(file="PCA_plot_gender_raw_chrXY.pdf")   ; p2; dev.off()
pdf(file="PCA_plot_gender_raw_chr1_22.pdf") ; p3; dev.off()

pdf(file="PCA_plot_gender_vst_all_chr.pdf") ; p4; dev.off()
pdf(file="PCA_plot_gender_vst_chrXY.pdf")   ; p5; dev.off()
pdf(file="PCA_plot_gender_vst_chr1_22.pdf") ; p6; dev.off()

pdf(file="PCA_plot_gender_rlog_all_chr.pdf") ; p7; dev.off()
pdf(file="PCA_plot_gender_rlog_chrXY.pdf")   ; p8; dev.off()
pdf(file="PCA_plot_gender_rlog_chr1_22.pdf") ; p9; dev.off()


# extract the legend and plot it separately
p0 <- plotPCA(dds_all_trans, intgroup=c("gender"))
leg <- get_legend(p0)
pdf(file="PCA_legend.pdf", width=3, height=3)
plot(leg)
dev.off()

# there's more code for adding sample names as labels on the plots, etc

################################################################################
################################################################################
################################################################################

# extract transformed values (Variance stabilizing transformation)
vst <- vst(dds, blind=FALSE)
rlog <- rlog(dds, blind=FALSE)
#head(assay(vst), 3)

ntd <- normTransform(dds)

library(vsn)

meanSdPlot(assay(ntd))    # Normalized counts transformation
meanSdPlot(assay(vst))    # Variance stabilizing transformation
meanSdPlot(assay(rlog))    # Regularized log transformation

# save other count matrices as well
write.csv(assay(ntd), file="normalized_counts_transformation_counts.csv")
write.csv(assay(vst), file="variance_stabilizing_transformation_counts.csv")
write.csv(assay(rlog), file="regularized_log_transformation_counts.csv")

#-----------------------------------------------------------------------------------------------#
# copied from 'run_Distribution_plot.R'

# raw counts (as they came from featureCounts())
data.raw <- counts_mat

# raw filtered counts (pre-filtering applied to keep only rows that have at least 10 reads total)
data.raw.filtered <- counts(dds)

# logarithm transformation => it will get rid of some extreme values. 
data.log2.on.raw <- log2(data.raw + 1)

# logarithm transformation => it will get rid of some extreme values. 
data.log2.on.filtered <- log2(data.raw.filtered + 1)

# variance-stabilizing transformation (VST), implemented in the DESeq package (Anders and Huber, 2010)
data.vst.on.raw <- vst(data.raw)

# variance-stabilizing transformation (VST), implemented in the DESeq package (Anders and Huber, 2010)
data.vst.on.filtered <- vst(data.raw.filtered)

# normalised counts (by DESeq)
data.norm <- counts_norm

# other DESeq2 transformations
# assay(ntd)
# assay(vst)
# assay(rlog)

# visualise 4 genes: 1st, 2nd, 14th and 18th; can't just take random genes as there are many with counts around 0
# data.raw            |     data.log2.on.raw        |     data.vst.on.raw         | dim= 58608    41
# data.raw.filtered   |     data.log2.on.filtered   |     data.vst.on.filtered    | dim= 52239    41

# data.norm     | dim= 52239    41

### FIRST PLOT
pdf(file="density_plots_DESeq2-p1.pdf", width=12, height=12)
par(mfrow=c(3,3))
plot(density(as.numeric(data.raw[1,])), main=paste0("raw - ", rownames(data.raw)[1]), cex.main=1,)
plot(density(as.numeric(data.log2.on.raw[1,])), main=paste0("log2 - ", rownames(data.log2.on.raw)[1]), cex.main=1)
plot(density(as.numeric(data.vst.on.raw[1,])), main=paste0("VST - ", rownames(data.vst.on.raw)[1]), cex.main=1)

plot(density(as.numeric(data.raw[14,])), main=paste0("raw - ", rownames(data.raw)[14]), cex.main=1,)
plot(density(as.numeric(data.log2.on.raw[14,])), main=paste0("log2 - ", rownames(log2.on.raw)[14]), cex.main=1)
plot(density(as.numeric(data.vst.on.raw[14,])), main=paste0("VST - ", rownames(data.vst.on.raw)[14]), cex.main=1)

plot(density(as.numeric(data.raw[18,])), main=paste0("raw - ", rownames(data.raw)[18]), cex.main=1,)
plot(density(as.numeric(data.log2.on.raw[18,])), main=paste0("log2 - ", rownames(log2.on.raw)[18]), cex.main=1)
plot(density(as.numeric(data.vst.on.raw[18,])), main=paste0("VST - ", rownames(data.raw)[18]), cex.main=1)
dev.off()

### SECOND PLOT
pdf(file="density_plots_DESeq2-p2.pdf", width=12, height=12)
par(mfrow=c(3,3))
plot(density(as.numeric(data.raw.filtered[1,])), main=paste0("raw (filtered) - ", rownames(data.raw.filtered)[1]), cex.main=1,)
plot(density(as.numeric(data.log2.on.filtered[1,])), main=paste0("log2 (filtered) - ", rownames(data.log2.on.filtered)[1]), cex.main=1)
plot(density(as.numeric(data.vst.on.filtered[1,])), main=paste0("VST (filtered) - ", rownames(data.vst.on.filtered)[1]), cex.main=1)

plot(density(as.numeric(data.raw[14,])), main=paste0("raw (filtered) - ", rownames(data.raw)[14]), cex.main=1,)
plot(density(as.numeric(data.log2.on.filtered[14,])), main=paste0("log2 (filtered) - ", rownames(data.log2.on.filtered)[14]), cex.main=1)
plot(density(as.numeric(data.vst.on.filtered[14,])), main=paste0("VST (filtered) - ", rownames(data.vst.on.filtered)[14]), cex.main=1)

plot(density(as.numeric(data.raw[18,])), main=paste0("raw (filtered) - ", rownames(data.raw)[18]), cex.main=1,)
plot(density(as.numeric(data.log2.on.filtered[18,])), main=paste0("log2 (filtered - ", rownames(data.log2.on.filtered)[18]), cex.main=1)
plot(density(as.numeric(data.vst.on.filtered[18,])), main=paste0("VST (filtered) - ", rownames(data.vst.on.filtered)[18]), cex.main=1)
dev.off()

### THIRD PLOT
pdf(file="density_plots_DESeq2-p3.pdf", width=12, height=12)
par(mfrow=c(3,3))
plot(density(as.numeric(data.norm[1,])), main=paste0("normalised (DESeq2) - ", rownames(data.norm)[1]), cex.main=1,)
plot(density(as.numeric(assay(vst)[1,])), main=paste0("transformed vst (DESeq2) - ", rownames(assay(vst))[1]), cex.main=1,)
plot(density(as.numeric(assay(rlog)[1,])), main=paste0("transformed R-log (DESeq2) - ", rownames(assay(rlog))[1]), cex.main=1,)

plot(density(as.numeric(data.norm[14,])), main=paste0("normalised (DESeq2) - ", rownames(data.norm)[14]), cex.main=1,)
plot(density(as.numeric(assay(vst)[14,])), main=paste0("transformed vst (DESeq2) - ", rownames(assay(vst))[14]), cex.main=1,)
plot(density(as.numeric(assay(rlog)[14,])), main=paste0("transformed R-log (DESeq2) - ", rownames(assay(rlog))[14]), cex.main=1,)

plot(density(as.numeric(data.norm[18,])), main=paste0("normalised (DESeq2) - ", rownames(data.norm)[18]), cex.main=1,)
plot(density(as.numeric(assay(vst)[18,])), main=paste0("transformed vst (DESeq2) - ", rownames(assay(vst))[18]), cex.main=1,)
plot(density(as.numeric(assay(rlog)[18,])), main=paste0("transformed Rilog (DESeq2) - ", rownames(assay(rlog))[18]), cex.main=1,)
dev.off()

cat("Finished!")
cat("Created: density_plots_DESeq2-p1.pdf")
cat("Created: density_plots_DESeq2-p2.pdf")
cat("Created: density_plots_DESeq2-p3.pdf")

#-----------------------------------------------------------------------------------------------#

