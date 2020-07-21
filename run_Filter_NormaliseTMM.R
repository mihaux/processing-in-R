# this script performs normalisation on raw counts data using trimmed mean of M values (TMM) [edgeR package]
# source: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
# source: https://seqqc.wordpress.com/2020/02/17/removing-low-count-genes-for-rna-seq-downstream-analysis/

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
#if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma"); library(limma)

# create a shortcut for the OneDrive directory where all files are stored
main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"      # on my mac
# main_dir <- "/Users/ummz/OneDrive - University of Leeds"                # on uni mac

#args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to .csv file with count data, 
       \n(2 - annotation) path to .csv annotation file and 
       \n(3 - output) path where output files should be stored", call.=FALSE)
}

# Example of usage: 
# Rscript run_Filter_NormaliseTMM.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/downstream/rerun_FINAL_July20/run_1/normalised_data/             

args <- c(paste0(main_dir, "/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_mod.csv"),
          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/run_1/normalised_data/"))

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load count data
df <- read.csv(args[1], row.names = 1, header = TRUE) # data.frame with counts only

# load annotation (clinical) data
anno <- read.csv(args[2], row.names = 1)

# add "ID_" to all rownames
rownames(anno) <- paste0("ID_", rownames(anno))

# create a vector with groups based on batch information
#group_temp <- gsub("one", 1, anno$Batch)
#group <- gsub("two", 2, group_temp)

# create a list-based data object called a DGEList
y <- DGEList(counts=df, group=anno$gender..1.male..2.female.)     # dim(y) => [1] 28278    41

# y$counts => a matrix counts containing the integer counts
# y$samples => a data.frame samples containing information about the samples or libraries

# The data.frame samples contains a column lib.size for the library size or sequencing depth for each sample.
# If not specified by the user, the library sizes will be computed from the column sums of the counts. 
# For classic edgeR the data.frame samples must also contain a column group, identifying the group membership of each sample.

### 1) Filtering

# As a rule of thumb, genes are dropped if they can’t possibly be expressed in all the samples for any of the conditions. 
# Users can set their own definition of genes being expressed. 
# Usually a gene is required to have a count of 5-10 in a library to be considered expressed in that library. 
# Users should also filter with count-per-million (CPM) rather than filtering on the counts directly, as the latter does not account for differences in library sizes between samples.

# filter out lowly expressed genes
# The filterByExpr() function keeps rows that have worthwhile counts in a minumum number of samples (15 samples in this case because the smallest group size is 15). 
# The function accesses the group factor contained in y in order to compute the minimum group size, but the filtering is performed independently of which sample belongs to which group so that no bias is introduced. 
# The group factor or the experimental design matrix can also be given directly to the filterByExpr function if not already set in the DGEList object. 
# It is also recommended to recalculate the library sizes of the DGEList object after the filtering, although the downstream analysis is robust to whether this is done or not.
keep <- filterByExpr(y)
y_filtered <- y[keep, , keep.lib.sizes=FALSE]       # dim(y_filtered) => [1] 16409    41

cat(dim(y)[1] - dim(y_filtered)[1], "genes were dropped. Filtered counts matrix dimension:", dim(y_filtered)[1], "x", dim(y_filtered)[2])

### 2) Normalization

# The calcNormFactors() function normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes. 
# The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples. 
# We call the product of the original library size and the scaling factor the effective library size. 
# The effective library size replaces the original library size in all downsteam analyses.

# TMM is recommended for most RNA-Seq data where the majority (more than half) of the genes are believed not differentially expressed between any pair of the samples. The following commands perform the TMM normalization and display the normalization factors.
y_normalised <- calcNormFactors(y_filtered)

# get normalized quantities for plotting etc.
# cpm() => Compute counts per million (CPM) or reads per kilobase per million (RPKM).
y_normalised_cpm <- cpm(y_normalised)

### SAVE RESULTS
info_1 <- tail(unlist(str_split(args[1], "/")), n=1)
info_2 <- unlist(str_split(info_1, "_"))[3:5]
info_3 <- paste(info_2[1], info_2[2], info_2[3], sep = "_")

# write filteres and normalised counts as .csv
write.csv(y_filtered$counts, file=paste0(args[3], "filtered_counts_", info_3, ".csv"))
write.csv(y_normalised_cpm, file=paste0(args[3], "normalised_cpm_counts_", info_3, ".csv"))

