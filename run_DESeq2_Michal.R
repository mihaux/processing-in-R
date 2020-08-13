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
library(ggplot2)

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

if (length(args)!=4) {
  stop("4 arguments must be supplied: 
       \n(1 - input) path to _x.csv file with count data, 
       \n(2 - annotation I) path to _x.csv annotation file I (clinical feature)
       \n(3 - annotation II) path to _x.csv annotation file II (slide scores)
       \nand (4 - output) path where output files should be stored", call.=FALSE)
}

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3

args <- c(paste0(main_dir, "/ANALYSES/run_12_Aug20/5_counting/PE_noXY/all_counts_dups_PE_noXY_mod.csv"),
          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/data/metadata/slide_scores/slide_scores_v6.csv"),
          paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/DESeq2_analysis/no_chrXY/"))

# (IN) /ANALYSES/run_12_Aug20/5_counting/PE_all/all_counts_dups_PE_all_mod.csv
# (OUT) /ANALYSES/run_12_Aug20/6_downstream/DESeq2_analysis/all_chr/

# (IN) /ANALYSES/run_12_Aug20/5_counting/PE_noXY/all_counts_dups_PE_noXY_mod.csv
# (OUT) /ANALYSES/run_12_Aug20/6_downstream/DESeq2_analysis/no_chrXY/

# Example of usage: 
# Rscript run_DESeq2.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_x.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED_x.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/downstream/rerun_FINAL_July20/run_1/DESeq2/            

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[4], sep="\n")
setwd(args[4])

# load count data (cols -> samples | rows -> genes)
df <- read.csv(args[1], row.names = 1, header = TRUE)     # data.frame with counts only

# NOTICE: It is important to never provide counts that were pre-normalized for sequencing depth/library size, 
# as the statistical model is most powerful when applied to un-normalized counts, 
# and is designed to account for library size differences internally.

# NOTICE: the counts table must be in the form of a matrix of integer values
# switch from data.frame to matrix and to integer
counts_mat <- as.matrix(df)                       # class: matrix | type: double  | dim: 28278    41
storage.mode(counts_mat) <- "integer"             # class: matrix | type: integer

# load annotation (clinical feature and slide scores) data which will be used for coldata
anno_1 <- read.csv(args[2], row.names = 1)
anno_2 <- read.csv(args[3], row.names = 1)

# NOTE: sample "14058" is missing in slide scores

# add "ID_" to all rownames
rownames(anno_1) <- paste0("ID_", rownames(anno_1))
rownames(anno_2) <- paste0("ID_", rownames(anno_2))

# make df colnames matchin the annotation rownames
colnames(counts_mat) <- str_replace(colnames(df), "X", "ID_")

# create DESeqDataSet object from Count matrix input
# NOTE: A DESeqDataSet object must have an associated design formula. 
# The design formula expresses the variables which will be used in modeling (and is used to estimate the dispersions and to estimate the log2 fold changes of the model.)

# NOTE: It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. 

# make sure the order is the same for counts_mat and anno and names correspond
all(rownames(anno_1) %in% colnames(counts_mat))
all(rownames(anno_1) == colnames(counts_mat))

all(rownames(anno_2) %in% colnames(counts_mat))
#all(rownames(anno_2) == colnames(counts_mat))  # doesn't work as the anno_2 is shorter due to lack of samples ID_14058

coldata_1 <- cbind(rownames(anno_1), anno_1[,c("visual.loss.at.BL..0.no..1.yes.", 
                                               "jaw.claudication.at.BL...0.no..1.yes.", 
                                               "ischaemic.features.at.BL...0.no..1.yes.", 
                                               "gender..1.male..2.female.", 
                                               "year.TAB.sample.was.collected", 
                                               "number.of.days.on.steroids.at.TAB", 
                                               "age.at.BL")])

names(coldata_1) <- c("sample", "visual_loss", "jaw_claudication", "ischaemic_features", "gender", "year_TAB_collected", "days_on_steroids", "age")

coldata_2 <- cbind(rownames(anno_2), anno_2[,c("GCA.present.",                          
                                               "Lymphocytic.infiltrate.in.media.",
                                               "Lymphocytic.infiltrate.in.intima.",      
                                               "Adventitia.pattern.",     
                                               "Media.pattern.",
                                               "Intima.pattern.",                        
                                               "Giant.cells.",                      
                                               "Infiltrate.around.vasa.vasorum.",                               
                                               "Media.destruction.",
                                               "Neoangiogenesis.",                    
                                               "Hyperplasia.",
                                               "Fibrosis.",                          
                                               "Oedema.",
                                               "Occlusion.grade.")])
  
names(coldata_2) <- c("sample", "GCA_present", "Lymphocytic_infiltrate_media", "Lymphocytic_infiltrate_intima", "Adventitia_pattern", "Media_pattern", "Intima_pattern", "Giant_cells", "Infiltrate_vasa", "Media_destruction", "Neoangiogenesis", "Hyperplasia", "Fibrosis", "Oedema", "Occlusion_grade")

# from slide scores (main)
# GCA.present   |   Giant.cells   |   Media.destruction   |   Occlusion.grade   |   Neoangiogenesis

# from slide scores (others)
# Lymphocytic.infiltrate.in.media | Lymphocytic.infiltrate.in.intima | Adventitia.pattern | Media.pattern
# Intima.pattern | Infiltrate.around.vasa.vasorum | Hyperplasia | Fibrosis | Oedema | Occlusion.grade 


# transform coldata feature to text (otherwise, there will be an error when running DESeq and also in PCA plots the group coloring will be treated as color scale)
coldata_1$gender <- str_replace(coldata_1$gender, pattern = "1", "male")
coldata_1$gender <-str_replace(coldata_1$gender, pattern = "2", "female")

coldata_1$visual_loss <- str_replace(coldata_1$visual_loss, pattern = "0", "no")
coldata_1$visual_loss <- str_replace(coldata_1$visual_loss, pattern = "1", "yes")

coldata_1$jaw_claudication <- str_replace(coldata_1$jaw_claudication, pattern = "0", "no")
coldata_1$jaw_claudication <- str_replace(coldata_1$jaw_claudication, pattern = "1", "yes")

coldata_1$ischaemic_features <- str_replace(coldata_1$ischaemic_features, pattern = "0", "no")
coldata_1$ischaemic_features <- str_replace(coldata_1$ischaemic_features, pattern = "1", "yes")


year_TAB_collected_new <- c()
#days_on_steroids_new <- c()
#age_new <- c()
for (i in 1:nrow(coldata_1)) {
  year_TAB_collected_new[i] <- paste0("_", coldata_1$year_TAB_collected[i], "_")
#  days_on_steroids_new[i] <- paste0("_", coldata_1$days_on_steroids[i], "_")
#  age_new[i] <- paste0("_", coldata_1$age[i], "_")
}

coldata_1$year_TAB_collected <- year_TAB_collected_new
#coldata_1$days_on_steroids <- days_on_steroids_new
#coldata_1$age <- age_new


# set up factors
coldata_1$gender <- factor(coldata_1$gender)
#coldata_1$visual_loss <- factor(coldata_1$visual_loss)

dds <- DESeqDataSetFromMatrix(countData = counts_mat, 
                              colData = coldata_1, 
                              design = ~ gender)

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

# set "female" as reference level (no reason for this)
dds$gender <- factor(dds$gender, levels = c("female","male"))

# perform DESeq analysis
dds <- DESeq(dds)
res <- results(dds)       # dim = 52239     6

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

# NOTE: PCA function needs a DESeqTransform class object => DESeqTransform(SummarizedExperiment)

# unnormalised
dds_all <- dds

# it takes too long, especially for rlog, so just load saved data (below)
# vst (variance stabilizing transformation)
vst_all     <- vst(dds, blind=FALSE)

# rlog (it takes quite long time to run)
rlog_all       <- rlog(dds, blind=FALSE)

# save transformed data and unnormalised as well
save(dds_all, file="Raw_DESeq_dataset_noXY.Rda")
save(vst_all, file="Normalised_DESeq_vst_dataset_noXY.Rda")
save(rlog_all, file="Normalised_DESeq_rlog_dataset_noXY.Rda")

