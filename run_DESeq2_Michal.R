# this script performs DESeq2 analysis on read counts data (output from featureCounts)
# created based DESeq2 vignette on https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# DESEq2 manual: https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf 

# The package DESeq2 provides methods to test for differential expression by use of negative binomial generalized linear models; 
# the estimates of dispersion and logarithmic fold changes incorporate data-driven prior distributions.

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr"); library(stringr)

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

#args <- c(paste0(main_dir, "/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_mod_x.csv"),
#          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
#          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/run_1/DESeq2/"))

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3

args <- c(paste0(main_dir, "/ANALYSES/comparison_with_Ian_results/rerun_5/featCounts/all_counts_dups_rr5_x.csv"),
          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/DESeq2/"))

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

coldata <- anno[,c("visual.loss.ever...0.no..1.yes." ,"gender..1.male..2.female.")]
names(coldata) <- c("visual_loss", "gender")

# transform coldata feature to text (not sure if it's needed)
coldata$visual_loss <- str_replace(coldata$visual_loss, pattern = "0", "no")
coldata$visual_loss <- str_replace(coldata$visual_loss, pattern = "1", "yes")

coldata$gender <- str_replace(coldata$gender, pattern = "1", "male")
coldata$gender <-str_replace(coldata$gender, pattern = "2", "female")

coldata$visual_loss <- factor(coldata$visual_loss)
coldata$gender <- factor(coldata$gender)

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~ gender)

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
dds <- dds[keep,]                     # 58608 - 52239 = 6369 rows dropped

# set "female" as reference level (no reason for this)
dds$gender <- factor(dds$gender, levels = c("female","male"))

# PCAplot to see if they cluster togther
# TODO


# perform DESeq analysis
dds <- DESeq(dds)
res <- results(dds)       # dim = 52239     6

# resultsNames(dds)  # [1] "Intercept"             "gender_male_vs_female"

mcols(res, use.names=TRUE)

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
plot.new()
plot.new()

plot(density(as.numeric(data.norm[14,])), main=paste0("normalised (DESeq2) - ", rownames(data.norm)[14]), cex.main=1,)
plot.new()
plot.new()

plot(density(as.numeric(data.norm[18,])), main=paste0("normalised (DESeq2) - ", rownames(data.norm)[18]), cex.main=1,)
plot.new()
plot.new()
dev.off()

cat("Finished!")
cat("Created: density_plots_DESeq2-p1.pdf")
cat("Created: density_plots_DESeq2-p2.pdf")
cat("Created: density_plots_DESeq2-p3.pdf")

#-----------------------------------------------------------------------------------------------#







