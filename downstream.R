# this script that to perform several analyses on read counts data (output from featureCounts)

# Summary of results to be generated:
# (1) removeBatchEffect() [edgeR package]
# (2) normalisation (trimmd mean of M method) [edgeR package]
# (3) hierarchical clustering
# (4) PCA
# (5) 

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db)
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma"); library(limma)

# INPUT (mode_I - many files: counts_merged.csv & stats_merged.csv | mode_II -> one file: all_counts_SE.csv & all_stats_SE.csv)
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_I
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_II

# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/paired-end/processed/mode_I
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/paired-end/processed/mode_II

# OUTPUT
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end
# /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/paired-end

#args <- commandArgs(trailingOnly = TRUE)

args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv",
          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/CiC_Clinical_data_FINAL.csv",
          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end")

cat("Example of usage: \n Rscript downstream.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/CiC_Clinical_data_FINAL.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end")

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

# if running for mode_II then, the colnames need to be changed
IDs <- sub(".Aligned.sortedByCoord.out.bam*", "", colnames(df))
IDs_final <- sub("X*", "", IDs)
colnames(df) <- IDs_final

# load annotation data
anno <- read.csv(args[2], row.names = 1)

# Visual loss at BL = reduced or lost vision (temporary or permanent) pre steroids or at BL assessment or AION or PION 
# Jaw claudication at BL = present pre steroids or at baseline assessment visit
# Ischaemic features at BL = jaw claudication, tongue claudication or visual loss (temporary or permanent) pre steroids or at baseline

# for columns 3-9: 0 - NO | 1 - YES
# for columns 5-6: 99 - not assessed
# for column 9: 2 - unknown
# for column 10: 1 - male | 2 - female

# select only necessary columns from anno
# [1] "Myriad" => not sure what it is ???                                    
# [2] "Batch"                                      
# [3] "visual_loss_at_BL"                  
# [4] "visual_loss_ever"                  
# [5] "AION_or_PION_BL"  
# [6] "AION_or_PION_ever"
# [7] "jaw_claudication_at_BL"            
# [8] "ischaemic_features_at_BL"          
# [9] "tongue_claudication"         
# [10] "gender"                        
# [11] "year_TAB_sample_was_collected"                    
# [12] "number_of_days_on_steroids_at_TAB"                
# [13] "number_of_days_between_TAB_and_BL_blood_sample" 

#---------------------------------------------#
### Principal Component Analysis (STANDARD) ###
#---------------------------------------------#
# perform scaling of data
dat <- scale(t(dat), scale=F)

#Therearetwobuilt-infunctionsinR toperformPCA:princompandprcomp.The former will perform eigen decomposition on the covariance matrix of the data, while the latter will perform singular value decomposition (SVD) on the data matrix. The latter is generally preferred. Eigen calculation on the covariance matrix (of size 22K by 22K in the original dimension of the data) can be time consuming. The calculation of SVD can be done more effectively on the data.


# you will have your own list here
symbols <- c("100287102", "653635", "102466751", "100302278")

# use mapIds method to obtain Entrez IDs
mapIds(org.Hs.eg.db, symbols, 'SYMBOL')




