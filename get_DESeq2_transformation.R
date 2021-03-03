# this script performs counts transformation (vst and rlog) from DESeq2 package 
# on raw read counts data (output from Salmon)

NOTE: this script has never been finished !!!


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); suppressMessages(library(DESeq2))

# get working directory to recognise the machine
w_dir <- getwd()

# create a shortcut for the OneDrive directory where all files are stored
if(startsWith(w_dir, "/Users/michal")){           
  main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"    # on my mac
} else if (startsWith(w_dir, "/Users/ummz")) {    
  main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds"      # on uni mac    
} else {
  print("Unrecognised machine.")
}

# define wheather output files should be saved or not [TRUE / FALSE]
output_save <- TRUE
          
# define directory with data (INPUT)
data_dir <- paste0(main_dir,"/ANALYSES_archived/run_13_Jan21/statistical_testing/input/")

# define directory for results (OUTPUT)
dir_out <- paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/statistical_testing/input_per_feature/")
setwd(dir_out)

# define directory with metadata (clinical and histological features)
dir_anno <- paste0(main_dir, "/data/metadata/outliers_excluded/")

# load raw counts
counts_raw_transcript <- read.csv(paste0(data_dir, "counts_raw_transcript-level.csv"), row.names = 1, header = TRUE)
counts_raw_gene <- read.csv(paste0(data_dir, "counts_raw_gene-level.csv"), row.names = 1, header = TRUE)

# load annotation files
anno_clinical <- read.csv(paste0(dir_anno, "cic_clinical_data_v2_summary_ORDERED_outliers_excluded.csv"), row.names = 1, header = TRUE)
anno_histological <- read.csv(paste0(dir_anno, "slide_scores_v6_outliers_excluded.csv"), row.names = 1, header = TRUE)

# add proper colnames to counts matrices
colnames(counts_raw_transcript) <- rownames(anno_clinical)
colnames(counts_raw_gene) <- rownames(anno_clinical)

# check if all col/row names are correct
#counts_raw_transcript[1:5,1:5]
#counts_raw_gene[1:5,1:5]
#anno_clinical[1:5,1:5]
#anno_histological[1:5,1:5]

# NOTICE: the counts table must be in the form of a matrix of integer values
# switch from data.frame to matrix and to integer
#counts_mat <- as.matrix(df)                       # class: matrix | type: double  | dim: 28278    41
#storage.mode(counts_mat) <- "integer"             # class: matrix | type: integer

counts_raw_transcript_fin <- as.matrix(counts_raw_transcript)
storage.mode(counts_raw_transcript_fin) <- "integer" 

counts_raw_gene_fin <- as.matrix(counts_raw_gene)
storage.mode(counts_raw_gene_fin) <- "integer" 


# it takes too long, especially for rlog, so just load saved data (below)
# vst (variance stabilizing transformation)
vst_transcript <- vst(counts_raw_transcript_fin, blind=FALSE)

# rlog (it takes quite long time to run)
rlog_all       <- rlog(dds, blind=FALSE)


# TODO:
# 1) load data and cread DESeq object
# 2) perfom the transformation (one per each feature)
# 3) save transformed matrices (one per each feature)
# 4) perform statistical testing [NOT HERE]
#
