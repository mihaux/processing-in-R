# FastQC postprocessing script that allows to generate several plots and tables to summarize the FastQC results 

# Summary of results to be generated:
# (1) plot with summary of all results for all QC metrics []
# (2) plot with GC content distribution (identical to the one from FastQC) for all samples in one  []
# (3) plot with percentages of each duplication level []
# (4) plots with total number of reads obtained for R1 and R2 []
# (5) .txt file with overepresented sequences []

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ngsReports", quietly = TRUE)) BiocManager::install("ngsReports"); library(ngsReports)
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2"); library(ggplot2)
if (!requireNamespace("dplyr", quietly = TRUE)) BiocManager::install("dplyr"); library(dplyr)
if (!requireNamespace("pander", quietly = TRUE)) BiocManager::install("pander"); library(pander)
if (!requireNamespace("xtable", quietly = TRUE)) BiocManager::install("xtable"); library(xtable)
if (!requireNamespace("gridExtra", quietly = TRUE)) BiocManager::install("gridExtra"); require(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
 
if (length(args)!=2) {
  stop("2 arguments must be supplied: (1 - input) path to directory with data and (2 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[1], sep="\n")

cat("Directory for results (OUT): ")
cat(args[2], sep="\n")

# select all zipped files and create an S4 object to store all results per sample
files_both <- list.files(args[1], pattern = "fastqc.zip$", full.names = TRUE)
fdl_both <- FastqcDataList(files_both)

# same as above but separately for R1 and R2
files_R1 <- list.files(args[1], pattern = "R1_001_fastqc.zip$", full.names = TRUE)
fdl_R1 <- FastqcDataList(files_R1)
files_R2 <- list.files(args[1], pattern = "R2_001_fastqc.zip$", full.names = TRUE)
fdl_R2 <- FastqcDataList(files_R2)

########## (1) plot with summary of all results for all QC metrics ##########




########## (2) plot with GC content distribution for all samples in one ##########


########## (3) plot with percentages of each duplication level ##########


########## (4) plots with total number of reads obtained for R1 and R2 ##########


########## (5) .txt file with overepresented sequences ##########

if(FALSE){
# create a table with the number of reads obtained for each sample (R1 and R2 separately)
# extract all read numbers
reads_R1 <- readTotals(fdl_R1)
reads_R2 <- readTotals(fdl_R2)

# extract all sample IDs (for R1 only as R2 are the same)
IDs <- substr(reads_R1$Filename, start = 1, stop = 5)
IDs <- gsub('_', "", IDs)

# change filenames to keep ID only (i.e. remove S.._L005 and 001) NOTICE: it is not needed as it can be changed later inside the plot object
for (x in 1:length(fdl_R1)) {
  file_name <- unique(fdl_R1[[x]]@Summary$Filename)
  file_name_new <- str_replace(file_name, ".fastq.gz", "_fastqc.zip")
  fdl_R1[[x]]@Summary$Filename <- gsub("_S\\d+_L005_R1_001", "_R1", fdl_R1[[x]]@Summary$Filename)
}

for (x in 1:length(fdl_R2)) {
  file_name <- unique(fdl_R2[[x]]@Summary$Filename)
  file_name_new <- str_replace(file_name, ".fastq.gz", "_fastqc.zip")
  fdl_R2[[x]]@Summary$Filename <- gsub("_S\\d+_L005_R2_001", "_R2", fdl_R2[[x]]@Summary$Filename)
}

# this can be done much simpler with:
new_leb_R1 <- c(rep("11026_R1_DV200", 11), rep("11028_R1", 11), rep("11037_R1_DV200", 11), rep("12223_R1_DV200", 11),
                rep("12330_R1", 11), rep("12331_R1", 11), rep("12388_R1", 11), rep("12389_R1", 11), rep("12849_R1_DV200", 11), 
                rep("12855_R1_DV200", 11), rep("12898_R1", 11), rep("12954_R1", 11), rep("14058_R1", 11), rep("14455_R1", 11), 
                rep("14912_R1", 11), rep("14914_R1_DV200", 11), rep("14915_R1_DV200", 11), rep("14916_R1", 11), rep("14919_R1", 11), 
                rep("14920_R1", 11), rep("14922_R1", 11), rep("14925_R1_DV200", 11), rep("14926_R1", 11), rep("14927_R1", 11), 
                rep("14929_R1", 11), rep("14930_R1", 11), rep("14931_R1_DV200", 11), rep("16082_R1", 11), rep("16098_R1", 11), 
                rep("16109_R1", 11), rep("16111_R1_DV200", 11), rep("16125_R1", 11), rep("16137_R1", 11), rep("16142_R1", 11), 
                rep("16153_R1", 11), rep("16154_R1", 11), rep("16167_R1", 11), rep("16187_R1_DV200", 11), rep("16648_R1_DV200", 11), 
                rep("16649_R1_DV200", 11), rep("8546_R1", 11))

new_leb_R2 <- c(rep("11026_R2_DV200", 11), rep("11028_R2", 11), rep("11037_R2_DV200", 11), rep("12223_R2_DV200", 11),
                rep("12330_R2", 11), rep("12331_R2", 11), rep("12388_R2", 11), rep("12389_R2", 11), rep("12849_R2_DV200", 11), 
                rep("12855_R2_DV200", 11), rep("12898_R2", 11), rep("12954_R2", 11), rep("14058_R2", 11), rep("14455_R2", 11), 
                rep("14912_R2", 11), rep("14914_R2_DV200", 11), rep("14915_R2_DV200", 11), rep("14916_R2", 11), rep("14919_R2", 11), 
                rep("14920_R2", 11), rep("14922_R2", 11), rep("14925_R2_DV200", 11), rep("14926_R2", 11), rep("14927_R2", 11), 
                rep("14929_R2", 11), rep("14930_R2", 11), rep("14931_R2_DV200", 11), rep("16082_R2", 11), rep("16098_R2", 11), 
                rep("16109_R2", 11), rep("16111_R2_DV200", 11), rep("16125_R2", 11), rep("16137_R2", 11), rep("16142_R2", 11), 
                rep("16153_R2", 11), rep("16154_R2", 11), rep("16167_R2", 11), rep("16187_R2_DV200", 11), rep("16648_R2_DV200", 11), 
                rep("16649_R2_DV200", 11), rep("8546_R2", 11))

# create a plot with PASS/WARN/FAIL Status of each module
sum_R1 <- plotSummary(fdl_R1); sum_R2 <- plotSummary(fdl_R2)
sum_R1$labels$x <- "QC metrics"; sum_R1$labels$y <- "Sample IDs"
sum_R2$labels$x <- "QC metrics"; sum_R2$labels$y <- "Sample IDs"

# add information about DV200 results
sum_R1$data$Filename <- new_leb_R1; sum_R2$data$Filename <- new_leb_R1

# TODO: need to figure out what are the best parameteres here, as by default the plot is too small, and loses quality when zoomed
#jpeg('summary_plot_R1.jpg', width = 880, height = 880, pointsize = 12, quality = 100, bg = "white",); sum_R1; dev.off()
#jpeg('summary_plot_R2.jpg'); sum_R2; dev.off()



# visualise total number of obtained Reads (unique and duplicated)
total_plot <- plotReadTotals(fdl_both)
total_plot$labels$x <- "Sample IDs"

to_be_mod <- as.vector(total_plot$data$Filename)
modified <- gsub("_S\\d+_L005", "", to_be_mod)
total_plot$data$Filename <- modified

total_plot <- total_plot + theme(text = element_text(size=7), 
                                 axis.title.x = element_text(size=9),
                                 axis.title.y = element_text(size=9))

duplic_plot_R1 <- plotDupLevels(fdl_R1)
duplic_plot_R2 <- plotDupLevels(fdl_R2)

# all GC plots in one figure 
# GC_plot_R1 <- plotGcContent(fdl_R1, plotType = "line",  gcType = "Transcriptome") + theme(legend.position="none")
# GC_plot_R2 <- plotGcContent(fdl_R2, plotType = "line",  gcType = "Transcriptome") + theme(legend.position="none")

# plot all GC plots together 
# legend:
# red: GC content per read
# blue: Theoretical distribution

# NOTICE: remove the legend (+ theme(legend.position="none")
# NOTICE: remove subtitles (p1$labels$subtitle = NULL)
# NOTICE: max 12 plots (grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, ncol = 3))

# WARN 
# r1_w_1 <- plotGcContent(fdl_R1$`11028_S4_L005_R1_001_fastqc.zip`, species = "Hsapiens", gcType = "Transcriptome") + theme(legend.position="none")
# r1_w_1$labels$subtitle = NULL
 
# jpeg('GC_content_R1_warn_I.jpg')
# grid.arrange(r1_w_1, r1_w_2, r1_w_3, r1_w_4, r1_w_5, r1_w_6, r1_w_7, r1_w_8, r1_w_9, r1_w_10, r1_w_11, r1_w_12, 
#               ncol = 3, top = textGrob("READ 1: WARN",gp=gpar(fontsize=14, font=2)))
# dev.off()

# FAIL 
# r1_f_1 <- plotGcContent(fdl_R1$`12331_S14_L005_R1_001_fastqc.zip`, species = "Hsapiens", gcType = "Transcriptome")
# r1_f_1$labels$subtitle = NULL


}

#################################
# save all plots and .txt files #
#################################






