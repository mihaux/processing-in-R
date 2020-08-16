# this script performs PCA analysis on read counts data (output from featureCounts)

# tuto on datacamp
# https://www.datacamp.com/community/tutorials/pca-analysis-r

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); suppressMessages(library(edgeR))
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma"); suppressMessages(library(limma))
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); suppressMessages(library(DESeq2))
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr"); library(stringr)
if (!requireNamespace("qqman", quietly = TRUE)) install.packages("qqman"); library(qqman)
library(dplyr)

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
# Rscript run_PCA_built-in.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/data/metadata/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/rerun_FINAL_July20/ANALYSES/downstream/rerun_FINAL_July20/run_1/PCA_plots             

#args <- c(paste0(main_dir, "/ANALYSES/rerun_FINAL_July20/run_1/featCounts_SE/all_counts_dups_run1_SE_mod.csv"),
#          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
#          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/run_1/PCA_plots"))

# /Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/comparison_with_Ian_results/rerun_5/featCounts/all_counts_nodups_rr5.csv

args <- c(paste0(main_dir, "/ANALYSES/comparison_with_Ian_results/rerun_5/featCounts/all_counts_dups_rr5_x.csv"),
          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/PCA_plots/"))

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load count data
df <- read.csv(args[1], row.names = 1, header = TRUE) # data.frame with counts only

# load annotation (clinical data)
anno <- read.csv(args[2], row.names = 1)

# load file with Transcripts by chromosome information
transcript_chr <- read.csv(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/table_TranscriptsByChromosome.csv"), header = TRUE, sep = ";", row.names = 1)

# add "ID_" to all rownames
rownames(anno) <- paste0("ID_", rownames(anno))

# make df colnames matchin the annotation rownames
colnames(df) <- str_replace(colnames(df), "X", "ID_")

# replace 1 - male and 2 - female in "gender..1.male..2.female." column
#anno[,15] <- str_replace(anno[,15], "1", "male")
#anno[,15] <- str_replace(anno[,15], "2", "female")

# NOTE: use the wide format in PCA, i.e. samples on rows and genes on columns (transpose if necessary)

# TODO: need to perform some normalisation or log-transformation and also removing 0 or rows (genes) which contain mostly zeros

# To perform (standard) principal component analysis (PCA), we use the wide format of the data. 
# We need to substract the column means from each column. 
# Effectively, the gene expressions are set to have zero mean for each probeset.

# perform scaling of data (and transpose)
df_scaled <- scale(t(df), scale=F)

# perform log-transformation
df_log2 <- log2(df + 1)

# perform some kind of normalisation (e.g. VST from DESeq2 package)
# switch from data.frame to matrix and to integer
counts_mat <- as.matrix(df)                       # class: matrix | type: double
storage.mode(counts_mat) <- "integer"             # class: matrix | type: integer

df_vst <- vst(counts_mat)

# There are two built-in functions in R to perform PCA:
# princomp() => it performs eigen decomposition on the covariance matrix of the data
# prcomp() => it performs singular value decomposition (SVD) on the data matrix. (this one is preferred) 

# Eigen calculation on a large covariance matrix can be time consuming. 
# The calculation of SVD can be done more effectively on the data.

# perform singular value decomposition (SVD) on the data matrix
res_scaled <- prcomp(df_scaled)
res_log2 <- prcomp(t(df_log2))
res_vst <- prcomp(t(df_vst))

# NOTE: by default, the function prcomp will center the columns of dat (we did it above for good practice). 

# res$sdev      => contains the standard deviation of the data explained by each principal component (from the largest to the smallest). 
# res$rotation  => a matrix of loadings, whose columns correspond to the loadings for each principal component (ordered from the first to the last, from the left). 
# res$center    => a vector of column means of the data (which is subtracted from each column of the data before the calculation of principal component is done).
# res$scale     =>
# res$x         => contains a matrix of principal components (columns), sometimes called scores.

# extract info about data from input filename
info_1 <- tail(unlist(str_split(args[1], "/")), n=1)
info_2 <- unlist(str_split(info_1, "_"))[3:5]
info_3 <- paste(info_2[1], info_2[2], info_2[3], sep = "_")

# investigate the results of the PCA
# make a scatter plot between the first and second loadings, and first and second PC’s
pdf(file=paste0(args[3], "scatter_plots_", info_3, ".pdf"), width=12, height=12)
par(mfrow=c(3,2), las=1)
plot(res_scaled$rotation[,1], res_scaled$rotation[,2], xlab = "Loading 1", ylab = "Loading 2", main="Loadings (scaled data)")
plot(res_scaled$x[,1], res_scaled$x[,2], xlab = "PC1", ylab = "PC2", main="Principal components (scaled data)")

plot(res_log2$rotation[,1], res_log2$rotation[,2], xlab = "Loading 1", ylab = "Loading 2", main="Loadings (log2 transformed data)")
plot(res_log2$x[,1], res_log2$x[,2], xlab = "PC1", ylab = "PC2", main="Principal components (log2 transformed data)")

plot(res_vst$rotation[,1], res_vst$rotation[,2], xlab = "Loading 1", ylab = "Loading 2", main="Loadings (vst normalised data)")
plot(res_vst$x[,1], res_vst$x[,2], xlab = "PC1", ylab = "PC2", main="Principal components (vst normalised data)")
dev.off()

cat("Created:", paste0(args[3], "scatter_plots_", info_3, ".pdf"))

# make another plot to see the variance explained by the different principal components
pdf(file=paste0(args[3], "variance_plots_", info_3, ".pdf"), width=12, height=12)
par(mfrow=c(3,2), las=1)
plot(res_scaled$sdev^2, type="h", xlab = "", ylab = "", main="Eigen values (scaled data)")
plot(cumsum(res_scaled$sdev^2)/sum(res_scaled$sdev^2), xlab = "Number of PC", ylab = "Cummulative proportion", main="Cummulative eigen values (scaled data)")

plot(res_log2$sdev^2, type="h", xlab = "", ylab = "", main="Eigen values (log2 transformed data)")
plot(cumsum(res_log2$sdev^2)/sum(res_log2$sdev^2), xlab = "Number of PC", ylab = "Cummulative proportion", main="Cummulative eigen values (log2 transformed data)")

plot(res_vst$sdev^2, type="h", xlab = "", ylab = "", main="Eigen values (vst normalised data)")
plot(cumsum(res_vst$sdev^2)/sum(res_vst$sdev^2), xlab = "Number of PC", ylab = "Cummulative proportion", main="Cummulative eigen values (vst normalised data)")
dev.off()

cat("Created:", paste0(args[3], "variance_plots_", info_3, ".pdf"))

# plot PC’s to see how they are related to some of the clinical phenotypes (male and female)
# pdf(file = "test.pdf", width=4, height=4)
pdf(file=paste0(args[3], "by_gender_plots_", info_3, ".pdf"), width=12, height=12)
par(mfrow=c(3,2), las=1)
plot(res_scaled$x[,1], res_scaled$x[,2], xlab = "PC1", ylab = "PC2", main="Labeled by gender (scaled data)", pch=16, cex=1.2, col=anno[,15])
plot(res_scaled$x[,3], res_scaled$x[,2], xlab = "PC3", ylab = "PC2", main="Labeled by gender (scaled data)", pch=16, cex=1.2, col=anno[,15])

plot(res_log2$x[,1], res_log2$x[,2], xlab = "PC1", ylab = "PC2", main="Labeled by gender (log2 transformed data)", pch=16, cex=1.2, col=anno[,15])
plot(res_log2$x[,3], res_log2$x[,2], xlab = "PC3", ylab = "PC2", main="Labeled by gender (log2 transformed data)", pch=16, cex=1.2, col=anno[,15])

plot(res_vst$x[,1], res_vst$x[,2], xlab = "PC1", ylab = "PC2", main="Labeled by gender (vst normalised data)", pch=16, cex=1.2, col=anno[,15])
plot(res_vst$x[,3], res_vst$x[,2], xlab = "PC3", ylab = "PC2", main="Labeled by gender (vst normalised data)", pch=16, cex=1.2, col=anno[,15])
dev.off()

cat("Created:", paste0(args[3], "by_gender_plots_", info_3, ".pdf"))

# plot legend only as separate plot
pdf(file = "legend.pdf", width=3, height=3)
par(mfrow=c(1,1), las=1)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", c("Male", "Female"), pch=16, pt.cex=2, cex=1.5, col = c(1,2))
mtext("Gender", cex=2)
dev.off()


cat("Created:", paste0(args[3], "by_legend_plot.pdf"))

cat("Finished!")

#### Sparse Principal Component Analysis ####
if(FALSE){

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

# PCA LOADINGS INVESTIGATION

# check loading in terms of their location 

# load positions list
all_positions <- read.csv("/Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/PCA_plots/all_positions_unique.csv")

# => how many on which chromosome and 
cat(length(which(transcript_chr$chromosomes == "chr1")), "on chr1")

col1 <- c()
col2 <- c()
for (i in 1:length(unique(transcript_chr$chromosomes))) {
  col1[i] <- as.character(unique(transcript_chr$chromosomes)[i])
  col2[i] <- length(which(transcript_chr$chromosomes == unique(transcript_chr$chromosomes)[i]))
}

nb_loadings <- data.frame(cbind(col1, col2))
rownames(nb_loadings) <- col1

# remove unecessary elements:
unnecessary_ele <- c("chr1_KI270713v1_random", "chr14_GL000194v1_random", "chr17_GL000205v2_random", 
                     "chr22_KI270734v1_random", "chr4_GL000008v2_random", "chrUn_GL000195v1",
                     "chrUn_GL000213v1", "chrUn_GL000218v1", "chrUn_GL000219v1", "chrUn_GL000224v1", 
                     "chrUn_KI270741v1")

nb_loadings_final <- nb_loadings[-which(rownames(nb_loadings) %in% unnecessary_ele),]

unnecessary_indices <- which(transcript_chr$chromosomes %in% unnecessary_ele)


# it's not in desired order, create a new one
nb_loadings_new <- rbind(nb_loadings_final[1,], nb_loadings_final[12,],nb_loadings_final[16:22,], 
      nb_loadings_final[2:11,], nb_loadings_final[13:15,], nb_loadings_final[23:25,])


pdf(file=paste0(args[3], "PCA_loadings_per_chromosome.pdf"), width=12, height=12)
barplot(as.numeric(as.vector(nb_loadings_new$col2)), main="Number of loadings per chromosome", 
        xlab="Chromosomes", names.arg=nb_loadings_new$col1, las=2)
dev.off()


# => plot them in chromosome order (chr1, chr2, chr3 etc.) => Manhattan plot for PC1, PC2, PC3 and PC4

# in my case: one plot per PC | make for PC1, PC2, PC3 and PC4
# SNP => list with ID (NM_)       => [1] "character"
# CHR => list with chromosomes    => [1] "integer"
# BP  => generate separate list for each chromosome   => [1] "integer"
# P   => PC1 loadnings            => [1] "numeric"

# generate BP
ranges <- as.numeric(as.vector(nb_loadings_new$col2))
#[1] 5827 4008 3274 2189 2643 3025 2703 2040 2255 2822 3443 2924 1105 1858 2060 2359 3242  874 3705 1668
#[21]  793 1258    2 2269  240
my_bp <- c(seq(1,ranges[1]), seq(1,ranges[2]), seq(1,ranges[3]), seq(1,ranges[4]), seq(1,ranges[5]),
           seq(1,ranges[6]), seq(1,ranges[7]), seq(1,ranges[8]), seq(1,ranges[9]), seq(1,ranges[10]),
           seq(1,ranges[11]), seq(1,ranges[12]), seq(1,ranges[13]), seq(1,ranges[14]), seq(1,ranges[15]),
           seq(1,ranges[16]), seq(1,ranges[17]), seq(1,ranges[18]), seq(1,ranges[19]), seq(1,ranges[20]),
           seq(1,ranges[21]), seq(1,ranges[22]), seq(1,ranges[23]), seq(1,ranges[24]), seq(1,ranges[25]))

my_snp <- as.character(as.vector(transcript_chr$IDs[-unnecessary_indices]))
my_chr <- as.vector(transcript_chr$chromosomes[-unnecessary_indices])
my_p1 <- as.vector(res_scaled$rotation[,1])
  
# load pc1_ordered 
pc1_ordered <- read.csv(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/PCA_plots/pc1_ordered.csv"))

to_be_matched <- data.frame(cbind(rownames(res_scaled$rotation), res_scaled$rotation[,1]))

# length(rownames(res_scaled$rotation)) => 58608 
# length(all_positions$x)               => 58586 (22 IDs fewer)

#diff_to_remove <- setdiff(rownames(res_scaled$rotation), all_positions$x)
#to_be_removed <- which(rownames(res_scaled$rotation) %in% diff_to_remove)
#my_p1_trimmed <- res_scaled$rotation[-to_be_removed,]
#write.csv(my_p1_trimmed, file = "pc1.csv")

# add ranking column to my_p1_trimmed[,1]





# length(rownames(my_p1_trimmed))   => 58586 OK



my_chr_fin_1 <- gsub("chr", "", my_chr)

# need to replace "chrM"  "chrX" and "chrY" with numerics: 23, 24, 25
my_chr_fin_2 <- gsub("M", "23", my_chr_fin_1)
my_chr_fin_3 <- gsub("X", "24", my_chr_fin_2)
my_chr_fin_4 <- gsub("Y", "25", my_chr_fin_3)
my_chr_final <- as.numeric(as.vector(my_chr_fin_4))


mygwas_list <- list(my_snp, as.integer(my_chr_final), as.integer(my_bp), as.numeric(my_p1[-unnecessary_indices])*(-1))
names(mygwas_list) <- c("SNP", "CHR", "BP", "P")

mygwas_list_bis <- list(as.character(as.vector(pc1_ordered$ID)), as.integer(my_chr_final), as.integer(my_bp), pc1_ordered$PC1*(-1))

#mygwas <- cbind(as.character(my_snp), as.vector(my_chr_final), as.integer(my_bp), as.numeric(my_p1[-unnecessary_indices]))
mygwas_fin <- as.data.frame(mygwas_list)
mygwas_fin_bis <- as.data.frame(mygwas_list_bis)

names(mygwas_fin) <- c("SNP", "CHR", "BP", "P")
names(mygwas_fin_bis) <- c("SNP", "CHR", "BP", "P")

manhattan(mygwas_fin, chr="CHR", bp="BP", snp="SNP", p="P", logp=FALSE, ylim = c(-0.05, 0.45))

manhattan(mygwas_fin_bis, chr="CHR", bp="BP", snp="SNP", p="P", logp=FALSE, ylim = c(-0.05, 0.45))

manhattan(mygwas_fin_bis, main = "Manhattan Plot - PC1", ylim = c(-0.05, 0.45), cex = 0.6, logp=FALSE,
          cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, 
          chrlabs = c(1:22, "M", "X", "Y"))

# => what exactly they are 







