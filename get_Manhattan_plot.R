# script that does Manhattan plot (written from scratch)
# it can be used for p-values but also for PCA loadings

# https://www.r-graph-gallery.com/101_Manhattan_plot.html
# https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R

# another method to create Manhattan plot on PCA loadings (TO BE MODIFIED)

# install (if necessary) and load package
library(DESeq2)
library(qqman)

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

setwd(paste0(main_dir, "/ANALYSES_archived/oct_nov/nov20_manhattan_plot"))


# run PCA using the built-in method and data with gene symbols
# define directory with data (INPUT)
data_dir_loadings <- paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/all_chr/INPUT_counts")

# load data RAW | VST | rlog (run for one data type at a time)
load(paste0(data_dir_loadings, "/Raw_DESeq_dataset_all.Rda")); dds <- dds_all
#load(paste0(data_dir_loadings,"/Normalised_DESeq_vst_dataset_all.Rda")); dds <- vst_all
#load(paste0(data_dir_loadings,"/Normalised_DESeq_rlog_dataset_all.Rda")); dds <- rlog_all

# run PCA using the built-in function
pca_dds <- prcomp(t(assay(dds)))

# retrieve PCA loadings
loadings_dds <- as.data.frame(pca_dds$rotation)

sink("rownames.txt")
for (i in 1:length(rownames(loadings_dds))) {
  cat(rownames(loadings_dds)[i], "\n")
}
sink()

# create a data.frame with
#'data.frame':	16470 obs. of  4 variables:
# $ SNP: chr  "rs1" "rs2" "rs3" "rs4" ...
# $ CHR: int  1 1 1 1 1 1 1 1 1 1 ...
# $ BP : int  1 2 3 4 5 6 7 8 9 10 ...
# $ P  : num  0.915 0.937 0.286 0.83 0.642 ...

# load sorted data
df_sorted <- read.csv("sorted.csv")

any(rownames(loadings_dds)[id_for_ordering] == df_sorted$Approved.symbol)
  
# create a data.frame with rownames, PC1 and PC2
df_final <- data.frame(genes=rownames(loadings_dds),
                       PC1=loadings_dds$PC1,
                       PC2=loadings_dds$PC2)

myResults_1 <- data.frame(SNP=rownames(loadings_dds)[id_for_ordering],
                          CHR=df_sorted$Chromosome,
                          BP=c(1:length(df_sorted$Approved.symbol)),
                          P=loadings_dds[id_for_ordering,]$PC1)

# change:
# => "unplaced"   to 23 
# => "X"          to 24
# => "X and Y"    to 25
# => "Y"          to 26

which(myResults_1$CHR == "unplaced")
which(myResults_1$CHR == "X")
which(myResults_1$CHR == "X and Y")
which(myResults_1$CHR == "Y")

# create a new vector for CHR
new_chr <- c(as.numeric(matrix(myResults_1$CHR[1:21741])), rep(23, 4), rep(24, 914), rep(25, 24), rep(26, 70))

# create a new vector with labels for chromosomes
chr_labels <- c(as.numeric(matrix(myResults_1$CHR[1:21741])), rep("unplaced", 4), rep("X", 914), rep("X and Y", 24), rep("Y", 70))


myResults_1_bis <- data.frame(SNP=rownames(loadings_dds)[id_for_ordering],
                              CHR=new_chr,
                              BP=c(1:length(df_sorted$Approved.symbol)),
                              P=as.numeric(loadings_dds[id_for_ordering,]$PC1))

myResults_2_bis <- data.frame(SNP=rownames(loadings_dds)[id_for_ordering],
                              CHR=new_chr,
                              BP=c(1:length(df_sorted$Approved.symbol)),
                              P=as.numeric(loadings_dds[id_for_ordering,]$PC2))

#manhattan(myResults_1_bis, chr="CHR", bp="BP", snp="SNP", p="P", logp=FALSE, chrlabs=as.character(chr_labels))
manhattan(myResults_1_bis, chr="CHR", bp="BP", snp="SNP", p="P", logp=FALSE, ylab="PC1 loadings", ylim=c(-0.01,0.8), main="Manhattan plot of PC1 loadings")
manhattan(myResults_2_bis, chr="CHR", bp="BP", snp="SNP", p="P", logp=FALSE, ylab="PC2 loadings", ylim=c(-0.8,0.4), main="Manhattan plot of PC2 loadings")

rankik <- c(1:length(df_sorted$Approved.symbol))

kot <- data.frame(ID=df_sorted$Approved.symbol, 
                  rank=rankik)

id_for_ordering <- match(kot$ID, rownames(loadings_dds))

head(rownames(loadings_dds)[id_for_ordering])

# get simple plot for PC1 and PC2 loadings
plot(myResults_1_bis$P, ylab="PC1 loadings", main="PC1 loadings - raw data")
plot(myResults_2_bis$P, ylab="PC2 loadings", main="PC2 loadings - raw data")



