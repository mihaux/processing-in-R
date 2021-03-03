# make plots of expression profile of a selected genes (using results from Salmon)
# data types: RAW and VST

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)
library(ggplot2)
library(plot.matrix)

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
dir_out <- paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/statistical_testing/expression_profiles/gene_level/")
setwd(dir_out)

# load data (dim = 31 896 x 40 - each)
df_raw <- read.csv(paste0(data_dir, "counts_raw_gene-level.csv"), row.names = 1)    
df_norm <- read.csv(paste0(data_dir, "counts_norm_gene-level.csv"), row.names = 1)
df_VST <- read.csv(paste0(data_dir, "counts_VST_gene-level.csv"), row.names = 1)

#############################################################################
# => check rows indices in df_raw (the row order is the same for all df_raw, df_norm and df_VST)

# SSTR1   |  ENSG00000139874
df_1 <- which(startsWith(rownames(df_raw), "ENSG00000139874"))    # 7835
rownames(df_raw)[df_1] # "ENSG00000139874"

# SSTR2   |  ENSG00000180616
df_2 <- which(startsWith(rownames(df_raw), "ENSG00000180616"))    # 14655
rownames(df_raw)[df_2]  # "ENSG00000180616"

# SSTR3   |  ENSG00000278195
df_3 <- which(startsWith(rownames(df_raw), "ENSG00000278195"))    # 31049
rownames(df_raw)[df_3]  # "ENSG00000278195"

# SSTR4   |  ENSG00000132671
df_4 <- which(startsWith(rownames(df_raw), "ENSG00000132671"))    # 6581
rownames(df_raw)[df_4] # "ENSG00000132671"

# SSTR5   |  ENSG00000162009
df_5 <- which(startsWith(rownames(df_raw), "ENSG00000162009"))    # 10578
rownames(df_raw)[df_5] # "ENSG00000162009

############################ _raw ############################
# get min, max, mean and median for each sample
df_raw_min <- apply(df_raw, 2, min)
df_raw_max <- apply(df_raw, 2, max)
df_raw_mean <- apply(df_raw, 2, mean)
df_raw_median <- apply(df_raw, 2, median)

# save results table
res_tab_1_raw <- data.frame(min=df_raw_min, max=df_raw_max, mean=df_raw_mean, median=df_raw_median)
if(output_save==TRUE){ write.csv(res_tab_1_raw, file = paste0("table_gene_raw.csv")) }

# create a matrix with all transcripts of interest
mat_1_raw <- rbind(df_raw[df_1,], df_raw[df_2,], df_raw[df_3,], df_raw[df_4,], df_raw[df_5,], 
                   df_raw_min, df_raw_max, df_raw_mean, df_raw_median)

rownames(mat_1_raw) <- c("SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5", "min", "max", "mean", "median")

# write a table 
write.csv(mat_1_raw, "genes_of_interest-counts_raw.csv")

# create a plot with all transcripts of interest and min, max, mean and median
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_gene_raw_all.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_raw), las=2, main="Expression profiles - raw - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_gene_raw_noMAX.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_raw[-7,]), las=2, main="Expression profiles - raw - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr and min.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_gene_raw_noMAXnoMIN.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_raw[-c(6,7),]), las=2, main="Expression profiles - raw - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

############################ _norm ############################
# get min, max, mean and median for each sample
df_norm_min <- apply(df_norm, 2, min)
df_norm_max <- apply(df_norm, 2, max)
df_norm_mean <- apply(df_norm, 2, mean)
df_norm_median <- apply(df_norm, 2, median)

# save results table
res_tab_1_norm <- data.frame(min=df_norm_min, max=df_norm_max, mean=df_norm_mean, median=df_norm_median)
if(output_save==TRUE){ write.csv(res_tab_1_norm, file = paste0("table_gene_norm.csv")) }

# create a matrix with all transcripts of interest
mat_1_norm <- rbind(df_norm[df_1,], df_norm[df_2,], df_norm[df_3,], df_norm[df_4,], df_norm[df_5,], 
                    df_norm_min, df_norm_max, df_norm_mean, df_norm_median)

rownames(mat_1_norm) <- c("SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5", "min", "max", "mean", "median")

# write a table 
write.csv(mat_1_norm, "genes_of_interest-counts_norm.csv")

# create a plot with all transcripts of interest and min, max, mean and median
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_gene_norm_all.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_norm), las=2, main="Expression profiles - norm - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_gene_norm_noMAX.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_norm[-7,]), las=2, main="Expression profiles - norm - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr and min.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_gene_norm_noMAXnoMIN.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_norm[-c(6,7),]), las=2, main="Expression profiles - norm - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

############################ _VST ############################
# get min, max, mean and median for each sample
df_VST_min <- apply(df_VST, 2, min)
df_VST_max <- apply(df_VST, 2, max)
df_VST_mean <- apply(df_VST, 2, mean)
df_VST_median <- apply(df_VST, 2, median)

# save results table
res_tab_1_VST <- data.frame(min=df_VST_min, max=df_VST_max, mean=df_VST_mean, median=df_VST_median)
if(output_save==TRUE){ write.csv(res_tab_1_VST, file = paste0("table_gene_VST.csv")) }

# create a matrix with all transcripts of interest
mat_1_VST <- rbind(df_VST[df_1,], df_VST[df_2,], df_VST[df_3,], df_VST[df_4,], df_VST[df_5,], 
                   df_VST_min, df_VST_max, df_VST_mean, df_VST_median)

rownames(mat_1_VST) <- c("SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5", "min", "max", "mean", "median")

# write a table 
write.csv(mat_1_VST, "genes_of_interest-counts_VST.csv")

# create a plot with all transcripts of interest and min, max, mean and median
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_gene_VST_all.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_VST), las=2, main="Expression profiles - VST - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_gene_VST_noMAX.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_VST[-7,]), las=2, main="Expression profiles - VST - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr and min.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_gene_VST_noMAXnoMIN.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_VST[-c(6,7),]), las=2, main="Expression profiles - VST - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }
