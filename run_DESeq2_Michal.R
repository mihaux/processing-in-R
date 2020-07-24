# this script performs DESeq2 analysis on read counts data (output from featureCounts)
# created based DESeq2 vignette on https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# DESEq2 manual: https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf 

# The package DESeq2 provides methods to test for differential expression by use of negative binomial generalized linear models; 
# the estimates of dispersion and logarithmic fold changes incorporate data-driven prior distributions.

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)

# if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
# if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr"); library(stringr)
# if (!requireNamespace("vsn", quietly = TRUE)) install.packages("vsn"); library(vsn)

# create a shortcut for the OneDrive directory where all files are stored
main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"      # on my mac
# main_dir <- "/Users/ummz/OneDrive - University of Leeds"                # on uni mac

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

#args <- c(paste0(main_dir, "/ANALYSES/comparison_with_Ian_results/rerun_5/featCounts/all_counts_dups_rr5.csv"),
#          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
 #         paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/DESeq2/"))


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


# load annotation (clinical) data
anno <- read.csv(args[2], row.names = 1)

# add "ID_" to all rownames
rownames(anno) <- paste0("ID_", rownames(anno))

# make df colnames matchin the annotation rownames
colnames(df) <- str_replace(colnames(df), "X", "ID_")



