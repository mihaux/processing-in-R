# this script performs PCA analysis on read counts data (output from featureCounts)

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma"); library(limma)

#args <- commandArgs(trailingOnly = TRUE)

args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/paired-end/normalised_cpm.csv",
          "/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary.csv",
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
df <- read.csv(args[1], row.names = 1, header = TRUE)

# if running for mode_II then, the colnames need to be changed
IDs <- sub(".Aligned.sortedByCoord.out.bam*", "", colnames(df))
IDs_final <- sub("X*", "", IDs)
colnames(df) <- IDs_final

# load annotation (clinical) data
anno <- read.csv(args[2], row.names = 1)

# TODO: need to perform some normalisation or log-transformation and also removing 0 or rows (genes) which contain mostly zeros
# NOTE: the data used is normalised

# To perform (standard) principal component analysis (PCA), we use the wide format of the data. 
# We need to substract the column means from each column. 
# Effectively, the gene expressions are set to have zero mean for each probeset.

# perform scaling of data
df_scaled <- scale(t(df), scale=F)

# There are two built-in functions in R to perform PCA:
# princomp() => it performs eigen decomposition on the covariance matrix of the data
# prcomp() => it performs singular value decomposition (SVD) on the data matrix. (this one is preferred) 

# Eigen calculation on a large covariance matrix can be time consuming. 
# The calculation of SVD can be done more effectively on the data.

# perform singular value decomposition (SVD) on the data matrix
res <- prcomp(df)

# NOTE: by default, the function prcomp will center the columns of dat (we did it above for good practice). 

# res$sdev      => contains the standard deviation of the data explained by each principal component (from the largest to the smallest). 
# res$rotation  => a matrix of loadings, whose columns corre- spond to the loadings for each principal component (ordered from the first to the last, from the left). 
# res$center    => a vector of column means of the data (which is subtracted from each column of the data before the calculation of principal component is done).
# res$scale     =>
# res$x         => contains a matrix of principal components (columns), sometimes called scores.

# investigate the results of the PCA
# make a scatter plot between the first and second loadings, and first and second PC’s
#par(mfrow=c(1,2), las=1)
#plot(res$rotation[,1], res$rotation[,2], xlab = "Loading 1", ylab = "Loading 2", main="Loadings")
#plot(res$x[,1], res$x[,2], xlab = "PC1", ylab = "PC2", main="Principal components")

# make another plot to see the variance explained by the different principal components
#par(mfrow=c(2,1), las=1)
#plot(res$sdev^2, type="h", xlab = "", ylab = "", main="Eigen values")
#plot(cumsum(res$sdev^2)/sum(res$sdev^2), xlab = "Number of PC", ylab = "Cummulative proportion", main="Cummulative eigen values")

# plot PC’s to see how they are related to some of the clinical phenotypes

#par(mfrow=c(1,2), las=1)

plot(res$x[,1], res$x[,2], xlab = "PC1", ylab = "PC2", main="Gender", col=anno[,15])
legend(50,40, c("ER-", "ER+"), col = c(1,2), pch=19)

plot(res$x[,3], res$x[,2], xlab = "PC3", ylab = "PC2", main="Gender", pch=19, cex=0.7, col=clinical[,8])
legend(50,40, c("ER-", "ER+"), col = c(1,2), pch=19)


if(FALSE){
#############################################
#### Sparse Principal Component Analysis ####
#############################################
# There are several packages available to perform sparse PCA 
# There is a function called arrayspc, however, we find that the results are not satisfactory (using default setting). 

install.packages("elasticnet")
library(elasticnet)

# run the built-in function spca() to get sparse loadings for the first four components, (NOTE: it may take a while!)
res2 <- spca(dat, K=4, rep(0.5,4))

spca.score1 <- dat%*%res2$loadings[,1]
spca.score2 <- dat%*%res2$loadings[,2]
spca.score3 <- dat%*%res2$loadings[,3]

# make similar plots as above to see the sparse loadings
plot(res2$loadings[,1], type="h", ylab="Loadings", main="Loadings for PC1")

# The proportion of variance explained for each principal component is contained in the component pev inside res2. 
# Note that there are just four sparse loadings that we calculated. 
# Had we calculated sparse loadings for a bigger number of loadings to calculate, then we can plot more points. 
# However, to calculate more loadings, it would take a very (!) long time due to the calculation of percentage of variance explained.

par(mfrow=c(2,1), las=1)
plot(res2$pev, type="h", ylab="Proportion", main="Variance Explained")
plot(cumsum(res2$pev), type="h", ylab = "Cummulative proportion", main="Variance explained", ylim=c(0,1))

# make the same plot to see the scores (from sparse loadings).
par(mfrow=c(1,2), las=1)
plot(spca.score1, spca.score2, xlab="SPC1", ylab="SPC2",
     main="ER Status (SPCA)", pch=19, cex=0.7, col=clinical[,8])
legend(50,40, c("ER-", "ER+"), col=c(1,2), pch=19)
plot(spca.score3, spca.score2, xlab="SPC3", ylab="SPC2",
     main="ER Status (SPCA)", pch=19, cex=0.7, col=clinical[,8])
legend(50,40, c("ER-", "ER+"), col=c(1,2), pch=19)

# OBSERVATION:
# the information on the scores are the same as the standard PCA. 
# However, the sparse PCA allows investigation on which genes/probesets contribute to the construction of the scores/PC.

}




