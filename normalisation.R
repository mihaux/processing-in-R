# this script performs normalisation on raw counts data using trimmed mean of M values (TMM) [edgeR package]
# source: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
#if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma"); library(limma)

# NOTICE: mode_I and mode_II results are identical, the only difference is in the order
# mode_II results are ordered in ascending way, but need to have modified colnames

# INPUT (mode_I - many files: counts_merged.csv & stats_merged.csv | mode_II -> one file: all_counts_SE.csv & all_stats_SE.csv)
#df_single <- read.csv("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv", row.names = 1)
#df_paired <- read.csv("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/paired-end/processed/mode_II/all_counts_PE.csv", row.names = 1)

# OUTPUT
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/paired-end

#args <- commandArgs(trailingOnly = TRUE)

args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv",
          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/CiC_Clinical_data_FINAL.csv",
          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end")

cat("Example of usage: \n Rscript normalisation.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/CiC_Clinical_data_FINAL.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end")

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to .csv file with count data, 
       \n(2 - annotation) path to .csv annotation file 
       \nand (3 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[1], sep="\n")

cat("Directory for results (OUT): ")
cat(args[3], sep="\n")

setwd(args[3])

# load count data
df <- read.csv(args[1], row.names = 1)

# the colnames need to be changed, as we are running for mode_II 
IDs <- sub(".Aligned.sortedByCoord.out.bam*", "", colnames(df))
IDs_final <- sub("X*", "", IDs)
colnames(df) <- IDs_final

# load annotation (clinical) data
anno <- read.csv(args[2], row.names = 1)

# create a vector with groups based on batch information
group_temp <- gsub("one", 1, anno$Batch)
group <- gsub("two", 2, group_temp)

# create a list-based data object called a DGEList
y <- DGEList(counts=df, group=group)
# dim(y) => [1] 28395    41

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

keep <- filterByExpr(y)
y_filtered <- y[keep, , keep.lib.sizes=FALSE]
# dim(y_filtered) => [1] 17830    41
# NOTICE: 28395 - 17830 = 10565 lowly expressed genes (rows) were dropped

# The filterByExpr() function keeps rows that have worthwhile counts in a minumum number of samples (15 samples in this case because the smallest group size is 15). 
# The function accesses the group factor contained in y in order to compute the minimum group size, but the filtering is performed independently of which sample belongs to which group so that no bias is introduced. 
# The group factor or the experimental design matrix can also be given directly to the filterByExpr function if not already set in the DGEList object. 
# It is also recommended to recalculate the library sizes of the DGEList object after the filtering, although the downstream analysis is robust to whether this is done or not.

### 2) Normalization

# The calcNormFactors() function normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes. 
# The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples. 
# We call the product of the original library size and the scaling factor the effective library size. 
# The effective library size replaces the original library size in all downsteam analyses.

# TMM is recommended for most RNA-Seq data where the majority (more than half) of the genes are believed not differentially expressed between any pair of the samples. The following commands perform the TMM normalization and display the normalization factors.
y_normalised <- calcNormFactors(y_filtered)
y_normalised$samples

# get normalized quantities for plotting etc.
# cpm() => Compute counts per million (CPM) or reads per kilobase per million (RPKM).
y_normalised_cpm <- cpm(y_normalised)



