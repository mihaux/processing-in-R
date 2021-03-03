# make plots of expression profile of a selected transcript (using results from Salmon)
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
dir_out <- paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/statistical_testing/expression_profiles/transcript_level/")
setwd(dir_out)

# load data (dim = 151 271 x 40 - each)
df_raw <- read.csv(paste0(data_dir, "counts_raw_transcript-level.csv"), row.names = 1)    
df_norm <- read.csv(paste0(data_dir, "counts_norm_transcript-level.csv"), row.names = 1)
df_VST <- read.csv(paste0(data_dir, "counts_VST_transcript-level.csv"), row.names = 1)

#----------------------------------------------------------------------------------------------------#
# load counts before filtering (in case, some of the transcripts of interest won't be present in df_raw/df_norm/df_VST, 
# then it would mean that they had to be filtered out beforehand, because thay had zeros across all the samples)
#df_no_filter <- read.csv("/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES_archived/run_13_Jan21/processed/counts_transcript-level_no_outliers.csv", row.names = 1)

# check in df_no_filter
#df_1_bis <- which(startsWith(rownames(df_no_filter), "ENST00000267377"))    # 86319
#rownames(df_no_filter)[df_1_bis] # "ENST00000267377.2"
#----------------------------------------------------------------------------------------------------#

#############################################################################
# => check rows indices in df_raw (the row order is the same for all df_raw, df_norm and df_VST)

# SSTR1-201   |  ENST00000267377.3
df_1 <- which(startsWith(rownames(df_raw), "ENST00000267377"))    # 76 200
rownames(df_raw)[df_1] # "ENST00000267377.2"

# SSTR2-201   |  ENST00000357585.4 
# SSTR2-202   |  ENST00000579323.5
df_2.1 <- which(startsWith(rownames(df_raw), "ENST00000357585"))  # 119 663
df_2.2 <- which(startsWith(rownames(df_raw), "ENST00000579323"))  # 119 662
rownames(df_raw)[df_2.1]  # "ENST00000357585.3"
rownames(df_raw)[df_2.2]  # "ENST00000579323.5"

# SSTR3-201   |  ENST00000610913.2
# SSTR3-202   |  ENST00000617123.1
df_3.1 <- which(startsWith(rownames(df_raw), "ENST00000610913"))  # 112 492
df_3.2 <- which(startsWith(rownames(df_raw), "ENST00000617123"))  # 112 493
rownames(df_raw)[df_3.1]  # "ENST00000610913.1"
rownames(df_raw)[df_3.2]  # "ENST00000617123.1"

# SSTR4-201   |  ENST00000255008.5
df_4 <- which(startsWith(rownames(df_raw), "ENST00000255008"))    # 127 948
rownames(df_raw)[df_4] # "ENST00000255008.4"

# SSTR5-201   |  ENST00000293897.5
df_5 <- which(startsWith(rownames(df_raw), "ENST00000293897"))    # 94 115
rownames(df_raw)[df_5] # "ENST00000293897.5"

############################ _raw ############################
# get min, max, mean and median for each sample
df_raw_min <- apply(df_raw, 2, min)
df_raw_max <- apply(df_raw, 2, max)
df_raw_mean <- apply(df_raw, 2, mean)
df_raw_median <- apply(df_raw, 2, median)

# save results table
res_tab_1_raw <- data.frame(min=df_raw_min, max=df_raw_max, mean=df_raw_mean, median=df_raw_median)
if(output_save==TRUE){ write.csv(res_tab_1_raw, file = paste0("table_transcripts_raw.csv")) }

# create a matrix with all transcripts of interest
mat_1_raw <- rbind(df_raw[df_1,], df_raw[df_2.1,], df_raw[df_2.2,], df_raw[df_3.1,], df_raw[df_3.2,], 
                   df_raw[df_4,], df_raw[df_5,], df_raw_min, df_raw_max, df_raw_mean, df_raw_median)

rownames(mat_1_raw) <- c("SSTR1-201", "SSTR2-201", "SSTR2-202", "SSTR3-201", "SSTR3-202", "SSTR4-201", "SSTR5-201", "min", "max", "mean", "median")

# write a table 
write.csv(mat_1_raw, "transcripts_of_interest-counts_raw.csv")

# create a plot with all transcripts of interest and min, max, mean and median
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_transcripts_raw_all.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_raw), las=2, main="Expression profiles - raw - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_transcripts_raw_noMAX.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_raw[-9,]), las=2, main="Expression profiles - raw - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr and min.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_transcripts_raw_noMAXnoMIN.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_raw[-c(8,9),]), las=2, main="Expression profiles - raw - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

############################ _norm ############################
# get min, max, mean and median for each sample
df_norm_min <- apply(df_norm, 2, min)
df_norm_max <- apply(df_norm, 2, max)
df_norm_mean <- apply(df_norm, 2, mean)
df_norm_median <- apply(df_norm, 2, median)

# save results table
res_tab_1_norm <- data.frame(min=df_norm_min, max=df_norm_max, mean=df_norm_mean, median=df_norm_median)
if(output_save==TRUE){ write.csv(res_tab_1_norm, file = paste0("table_transcripts_norm.csv")) }

# create a matrix with all transcripts of interest
mat_1_norm <- rbind(df_norm[df_1,], df_norm[df_2.1,], df_norm[df_2.2,], df_norm[df_3.1,], df_norm[df_3.2,], 
                    df_norm[df_4,], df_norm[df_5,], df_norm_min, df_norm_max, df_norm_mean, df_norm_median)

rownames(mat_1_norm) <- c("SSTR1-201", "SSTR2-201", "SSTR2-202", "SSTR3-201", "SSTR3-202", "SSTR4-201", "SSTR5-201", "min", "max", "mean", "median")

# write a table 
write.csv(mat_1_norm, "transcripts_of_interest-counts_norm.csv")

# create a plot with all transcripts of interest and min, max, mean and median
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_transcripts_norm_all.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_norm), las=2, main="Expression profiles - norm - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_transcripts_norm_noMAX.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_norm[-9,]), las=2, main="Expression profiles - norm - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr and min.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_transcripts_norm_noMAXnoMIN.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_norm[-c(8,9),]), las=2, main="Expression profiles - norm - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

############################ _VST ############################
# get min, max, mean and median for each sample
df_VST_min <- apply(df_VST, 2, min)
df_VST_max <- apply(df_VST, 2, max)
df_VST_mean <- apply(df_VST, 2, mean)
df_VST_median <- apply(df_VST, 2, median)

# save results table
res_tab_1_VST <- data.frame(min=df_VST_min, max=df_VST_max, mean=df_VST_mean, median=df_VST_median)
if(output_save==TRUE){ write.csv(res_tab_1_VST, file = paste0("table_transcripts_VST.csv")) }

# create a matrix with all transcripts of interest
mat_1_VST <- rbind(df_VST[df_1,], df_VST[df_2.1,], df_VST[df_2.2,], df_VST[df_3.1,], df_VST[df_3.2,], 
                   df_VST[df_4,], df_VST[df_5,], df_VST_min, df_VST_max, df_VST_mean, df_VST_median)

rownames(mat_1_VST) <- c("SSTR1-201", "SSTR2-201", "SSTR2-202", "SSTR3-201", "SSTR3-202", "SSTR4-201", "SSTR5-201", "min", "max", "mean", "median")

# write a table 
write.csv(mat_1_VST, "transcripts_of_interest-counts_VST.csv")

# create a plot with all transcripts of interest and min, max, mean and median
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_transcripts_VST_all.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_VST), las=2, main="Expression profiles - VST - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_transcripts_VST_noMAX.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_VST[-9,]), las=2, main="Expression profiles - VST - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr and min.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("plot_SSTR1-5_transcripts_VST_noMAXnoMIN.png")) }
par(mar=c(6, 5, 6, 5))
plot(t(mat_1_VST[-c(8,9),]), las=2, main="Expression profiles - VST - all", ylab="", xlab="")
if(output_save==TRUE){ dev.off() }
