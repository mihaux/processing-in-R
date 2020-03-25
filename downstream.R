# this script that to perform several analyses on read counts data (output from featureCounts)

# Summary of results to be generated:
# (1) removeBatchEffect() [edgeR package]
# (2) normalisation (trimmd mean of M method) [edgeR package]
# (3) hierarchical clustering
# (4) PCA
# (5) 

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)

# INPUT (mode_I - many files | mode_II -> one file)
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_I
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_II

# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/paired-end/processed/mode_I
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/paired-end/processed/mode_II

# OUTPUT
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/paired-end

#args <- commandArgs(trailingOnly = TRUE)

args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_II", 
          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end")

cat("Example of usage: \n Rscript downstream.R /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_IV/1_quality_control/report /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_IV/1_quality_control/postprocessed")

if (length(args)!=2) {
  stop("2 arguments must be supplied: \n(1 - input) path to directory with data and \n(2 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[1], sep="\n")

cat("Directory for results (OUT): ")
cat(args[2], sep="\n")

setwd(args[2])

# select all zipped files and create an S4 object to store all results per sample
files_both <- list.files(args[1], pattern = "fastqc.zip$", full.names = TRUE)
fdl_both <- FastqcDataList(files_both)

# same as above but separately for R1 and R2
files_R1 <- list.files(args[1], pattern = "_R1.fastqc.zip$", full.names = TRUE)
fdl_R1 <- FastqcDataList(files_R1)
files_R2 <- list.files(args[1], pattern = "_R2.fastqc.zip$", full.names = TRUE)
fdl_R2 <- FastqcDataList(files_R2)