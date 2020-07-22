# this script performs DESeq2 analysis on read counts data (output from featureCounts) [adapted mostly from Ian's script]
# use 'run_PCA_DESeq2.R' for running PCA analysis afterwards

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr"); library(stringr)
if (!requireNamespace("vsn", quietly = TRUE)) install.packages("vsn"); library(vsn)

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

args <- c(paste0(main_dir, "/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_mod.csv"),
          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/run_1/DESeq2/"))

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3

# Example of usage: 
# Rscript run_PCA_DESeq2.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/downstream/rerun_FINAL_July20/run_1/DESeq2            

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load count data
df <- read.csv(args[1], row.names = 1, header = TRUE)     # data.frame with counts only

# load annotation (clinical) data
anno <- read.csv(args[2], row.names = 1)

# switch from data.frame to matrix and to integer
counts_mat <- as.matrix(df)                       # class: matrix | type: double  | dim: 28278    41
storage.mode(counts_mat) <- "integer"             # class: matrix | type: integer

# filter out genes that ... (adapted from Ian, not sure how it is filtered)
counts_mat_filtered <- counts_mat[rowSums(counts_mat > 3) > 2,]   # | dim: 21577    41  (6701 genes filtered out)

cat(dim(counts_mat)[1] - dim(counts_mat_filtered)[1], "genes were dropped. Filtered counts matrix dimension:", dim(counts_mat_filtered)[1], "x", dim(counts_mat_filtered)[2])

# define conditions for experiment ( gender as conditions, replace 1-"male" and 2-"female")
sampleConditions_temp <- str_replace(anno[,15], "1", "male")    
sampleConditions <- str_replace(sampleConditions_temp, "2", "female")

# Sample set name
conditionGroupsForDataFirstIsReference = c("male","female")

# Create experiment design model
sampleTable <- data.frame(sample=colnames(df), condition=sampleConditions)

# Set row names in model design to sample names
rownames(sampleTable) <- sampleTable$sample

# Check all is well with names
all(rownames(sampleTable) %in% colnames(counts_mat_filtered))

# put data from counts file in to correct format and set factors and levels
dds <- DESeqDataSetFromMatrix(countData=counts_mat_filtered, colData=sampleTable, design=~condition)

colData(dds)$condition <- factor(colData(dds)$condition, levels=conditionGroupsForDataFirstIsReference)

# set up the dds object
featureData <- data.frame(gene=rownames(counts_mat_filtered))
head(featureData)

mcols(dds)<- DataFrame(mcols(dds), featureData)
mcols(dds)

#does the work of the differential expression analysis
dds <- DESeq(dds)

resultsNames(dds)

# make a number of result objects based on their different analysis methods
res <- results(dds)

summary(res)

### SAVE RESULTS
info_1 <- tail(unlist(str_split(args[1], "/")), n=1)
info_2 <- unlist(str_split(info_1, "_"))[3:5]
info_3 <- paste(info_2[1], info_2[2], info_2[3], sep = "_")

# write results names, summary and data
write.csv(resultsNames(dds), file=paste0(args[3], "Analysis_names_", info_3, ".csv"))
sink(file=paste0(args[3], "Res_summary_", info_3, ".txt")); summary(res); sink()
save(dds, file=paste0(args[3], "DESeqDataSet_", info_3, ".Rda"))
write.csv(ras.data.frame(res), file=paste0(args[3], "standard_analysis_Deseq2_", info_3, ".csv"))

dev.copy(png, paste0(args[3], "standard_analysis_MA_Plot_Deseq2_", info_3, ".png"))
plotMA(res, ylim=c(-2,2), main='STD')
dev.off()

# get only significant results (with pvalue < 0.01)
resSig <- subset(res, padj < 0.01)

write.csv(as.data.frame(resSig), file=paste0(args[3], "standard_analysis_sig_Deseq2_", info_3, ".csv"))
save(res, file=paste0(args[3], "DESeq_standard_Result_dataset_", info_3, ".Rda"))
save(resSig, file=paste0(args[3], "DESeq_standard_Result_dataset_sig_", info_3, ".Rda"))

# normalise counts
nt <- normTransform(dds)
save(nt, file=paste0(args[3], "Normalised_DESeq_normal_transformation_dataset_", info_3, ".Rda"))

dev.copy(png, paste0(args[3], "Normalised_DESeq_normal_transformation_mean_Sd_Plot_", info_3, ".png"))
meanSdPlot(assay(nt))
dev.off()

head(assay(nt), 10)




