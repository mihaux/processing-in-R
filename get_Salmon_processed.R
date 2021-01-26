# make plots of expression profile of a selected gene / transcript
# data types: Raw | Normalised_rlog | Normalised_vst

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); library(DESeq2)
library(ggplot2)

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
output_save <- FALSE

# define directory with data (INPUT)
data_dir_1 <- paste0(main_dir,"/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/FINAL")

# single-end:
#data_dir <- paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/SE/DESeq2_analysis/all_chr/INPUT_counts")

# paired-end:
data_dir_2 <- paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/all_chr/INPUT_counts")

# define directory for results (OUTPUT)
dir_out <- paste0(main_dir, "/ANALYSES/nov20_SSTR2-expr-profile/vst/")
setwd(dir_out)

# load data RAW | VST | rlog (run for one data type at a time)

# NOTE: the use of raw data does not make much sense here
#load(paste0(data_dir_1, "/Raw_DESeq_dataset_all.Rda")); dds_1 <- dds_all # 52 239 rows, RefSeq annotation
#load(paste0(data_dir_2, "/Raw_DESeq_dataset_all.Rda")); dds_2 <- dds_all # 26 486 rows, GeneSymbol annotation

#load(paste0(data_dir_1,"/Normalised_DESeq_vst_dataset_all.Rda")); dds_1 <- vst_all
#load(paste0(data_dir_2,"/Normalised_DESeq_vst_dataset_all.Rda")); dds_2 <- vst_all

load(paste0(data_dir_1,"/Normalised_DESeq_rlog_dataset_all.Rda")); dds_1 <- rlog_all
load(paste0(data_dir_2,"/Normalised_DESeq_rlog_dataset_all.Rda")); dds_2 <- rlog_all

# define running ID (either "raw", "vst" pr "rlog")
run_id <- "rlog"

#############################################################################
# => in dds_1, look for NM_001049, NM_001050, NM_001051, NM_001052, NM_001053
#############################################################################

# get min, max, mean and median for each sample
dds_1_min <- apply(assay(dds_1), 2, min)
dds_1_max <- apply(assay(dds_1), 2, max)
dds_1_mean <- apply(assay(dds_1), 2, mean)
dds_1_median <- apply(assay(dds_1), 2, median)

# save results table
res_tab_1 <- data.frame(min=dds_1_min, max=dds_1_max, mean=dds_1_mean, median=dds_1_median)
if(output_save==TRUE){ write.csv(res_tab_1, file = paste0("table_transcripts_", run_id,".csv")) }

dds_1_1 <- which(startsWith(rownames(dds_1), "NM_001049")) # 37043
dds_1_2 <- which(startsWith(rownames(dds_1), "NM_001050")) # 44359
dds_1_3 <- which(startsWith(rownames(dds_1), "NM_001051")) # 51880
dds_1_4 <- which(startsWith(rownames(dds_1), "NM_001052")) # 49679
dds_1_5 <- which(startsWith(rownames(dds_1), "NM_001053")) # 40277

# create a matrix with all genes of interest
mat_1 <- cbind(assay(dds_1)[dds_1_1,], 
               assay(dds_1)[dds_1_2,], 
               assay(dds_1)[dds_1_3,], 
               assay(dds_1)[dds_1_4,], 
               assay(dds_1)[dds_1_5,])

colnames(mat_1) <- c("NM_001049", "NM_001050", "NM_001051", "NM_001052", "NM_001053")
  
# create a heatmap with all 5 genes of interest and min, max, mean and median
if(output_save==TRUE){ png(file = paste0("heatmap_SSTR1-5_transcripts_", run_id, "_all.png")) }
heatmap(cbind(mat_1, dds_1_min, dds_1_max, dds_1_mean, dds_1_median), Rowv = NA, Colv = NA, cexRow=0.8, cexCol=0.8)
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("heatmap_SSTR1-5_transcripts_", run_id, "_noMAX.png")) }
heatmap(cbind(mat_1, dds_1_min, dds_1_mean, dds_1_median), Rowv = NA, Colv = NA, cexRow=0.8, cexCol=0.8)
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr and min.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("heatmap_SSTR1-5_transcripts_", run_id, "_noMAXnoMIN.png")) }
heatmap(cbind(mat_1, dds_1_mean, dds_1_median), Rowv = NA, Colv = NA, cexRow=0.8, cexCol=0.8)
if(output_save==TRUE){ dev.off() }

# make histograms of gene expression values per sample
png(paste0("histograms_transcripts_per_sample_", run_id, "_part-I.png"))
hist(assay(dds_1)[,1], xlab = colnames(assay(dds_1))[1], main = "Distribution of gene expression per sample")
dev.off()

png(paste0("histograms_transcripts_per_sample_", run_id, "_part-II.png"))
par(mfrow=c(4,5))
hist(assay(dds_1)[,2], xlab = colnames(assay(dds_1))[2], main = NULL)
hist(assay(dds_1)[,3], xlab = colnames(assay(dds_1))[3], main = NULL)
hist(assay(dds_1)[,4], xlab = colnames(assay(dds_1))[4], main = NULL)
hist(assay(dds_1)[,5], xlab = colnames(assay(dds_1))[5], main = NULL)
hist(assay(dds_1)[,6], xlab = colnames(assay(dds_1))[6], main = NULL)
hist(assay(dds_1)[,7], xlab = colnames(assay(dds_1))[7], main = NULL)
hist(assay(dds_1)[,8], xlab = colnames(assay(dds_1))[8], main = NULL)
hist(assay(dds_1)[,9], xlab = colnames(assay(dds_1))[9], main = NULL)
hist(assay(dds_1)[,10], xlab = colnames(assay(dds_1))[10], main = NULL)
hist(assay(dds_1)[,11], xlab = colnames(assay(dds_1))[11], main = NULL)
hist(assay(dds_1)[,12], xlab = colnames(assay(dds_1))[12], main = NULL)
hist(assay(dds_1)[,13], xlab = colnames(assay(dds_1))[13], main = NULL)
hist(assay(dds_1)[,14], xlab = colnames(assay(dds_1))[14], main = NULL)
hist(assay(dds_1)[,15], xlab = colnames(assay(dds_1))[15], main = NULL)
hist(assay(dds_1)[,16], xlab = colnames(assay(dds_1))[16], main = NULL)
hist(assay(dds_1)[,17], xlab = colnames(assay(dds_1))[17], main = NULL)
hist(assay(dds_1)[,18], xlab = colnames(assay(dds_1))[18], main = NULL)
hist(assay(dds_1)[,19], xlab = colnames(assay(dds_1))[19], main = NULL)
hist(assay(dds_1)[,20], xlab = colnames(assay(dds_1))[20], main = NULL)
hist(assay(dds_1)[,21], xlab = colnames(assay(dds_1))[21], main = NULL)
dev.off()

png(paste0("histograms_transcripts_per_sample_", run_id, "_part-III.png"))
par(mfrow=c(4,5))
hist(assay(dds_1)[,22], xlab = colnames(assay(dds_1))[22], main = NULL)
hist(assay(dds_1)[,23], xlab = colnames(assay(dds_1))[23], main = NULL)
hist(assay(dds_1)[,24], xlab = colnames(assay(dds_1))[24], main = NULL)
hist(assay(dds_1)[,25], xlab = colnames(assay(dds_1))[25], main = NULL)
hist(assay(dds_1)[,26], xlab = colnames(assay(dds_1))[26], main = NULL)
hist(assay(dds_1)[,27], xlab = colnames(assay(dds_1))[27], main = NULL)
hist(assay(dds_1)[,28], xlab = colnames(assay(dds_1))[28], main = NULL)
hist(assay(dds_1)[,29], xlab = colnames(assay(dds_1))[29], main = NULL)
hist(assay(dds_1)[,30], xlab = colnames(assay(dds_1))[30], main = NULL)
hist(assay(dds_1)[,31], xlab = colnames(assay(dds_1))[31], main = NULL)
hist(assay(dds_1)[,32], xlab = colnames(assay(dds_1))[32], main = NULL)
hist(assay(dds_1)[,33], xlab = colnames(assay(dds_1))[33], main = NULL)
hist(assay(dds_1)[,34], xlab = colnames(assay(dds_1))[34], main = NULL)
hist(assay(dds_1)[,35], xlab = colnames(assay(dds_1))[35], main = NULL)
hist(assay(dds_1)[,36], xlab = colnames(assay(dds_1))[36], main = NULL)
hist(assay(dds_1)[,37], xlab = colnames(assay(dds_1))[37], main = NULL)
hist(assay(dds_1)[,38], xlab = colnames(assay(dds_1))[38], main = NULL)
hist(assay(dds_1)[,39], xlab = colnames(assay(dds_1))[39], main = NULL)
hist(assay(dds_1)[,40], xlab = colnames(assay(dds_1))[40], main = NULL)
hist(assay(dds_1)[,41], xlab = colnames(assay(dds_1))[41], main = NULL)
dev.off()

#############################################################################
# => in dds_2, look for SSTR1, SSTR2, SSTR3, SSTR4, SSTR5
#############################################################################

# get min, max, mean and median for each sample
dds_2_min <- apply(assay(dds_2), 2, min)
dds_2_max <- apply(assay(dds_2), 2, max)
dds_2_mean <- apply(assay(dds_2), 2, mean)
dds_2_median <- apply(assay(dds_2), 2, median)

# save results table
res_tab_2 <- data.frame(min=dds_2_min, max=dds_2_max, mean=dds_2_mean, median=dds_2_median)
if(output_save==TRUE){ write.csv(res_tab_2, file = paste0("table_genes_", run_id,".csv")) }

dds_2_1 <- which(startsWith(rownames(dds_2), "SSTR1")) # 7321
dds_2_2 <- which(startsWith(rownames(dds_2), "SSTR2")) # 11356
dds_2_3 <- which(startsWith(rownames(dds_2), "SSTR3")) # 16787
dds_2_4 <- which(startsWith(rownames(dds_2), "SSTR4")) # 15559
dds_2_5 <- which(startsWith(rownames(dds_2), "SSTR5")) # 9008 and 9009

# create a matrix with all genes of interest
mat_2 <- cbind(assay(dds_2)[dds_2_1,], 
               assay(dds_2)[dds_2_2,], 
               assay(dds_2)[dds_2_3,], 
               assay(dds_2)[dds_2_4,], 
               #assay(dds_2)[dds_2_5[1],],  # => SSTR5 antisense RNA 1
               assay(dds_2)[dds_2_5[2],])

colnames(mat_2) <- c("SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5")

# create a heatmap with all 5 genes of interest and min, max, mean and median
if(output_save==TRUE){ png(file = paste0("heatmap_SSTR1-5_genes_", run_id, "_all.png")) }
heatmap(cbind(mat_2, dds_2_min, dds_2_max, dds_2_mean, dds_2_median), Rowv = NA, Colv = NA, cexRow=0.8, cexCol=0.8)
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("heatmap_SSTR1-5_genes_", run_id, "_noMAX.png")) }
heatmap(cbind(mat_2, dds_2_min, dds_2_mean, dds_2_median), Rowv = NA, Colv = NA, cexRow=0.8, cexCol=0.8)
if(output_save==TRUE){ dev.off() }

# same as above, but without max.expr and min.expr as it's too extreme
if(output_save==TRUE){ png(file = paste0("heatmap_SSTR1-5_genes_", run_id, "_noMAXnoMIN.png")) }
heatmap(cbind(mat_2, dds_2_mean, dds_2_median), Rowv = NA, Colv = NA, cexRow=0.8, cexCol=0.8)
if(output_save==TRUE){ dev.off() }

par(mfrow=c(1,1))
hist(assay(dds_2)[,1], xlab = colnames(assay(dds_2))[1], main = NULL)

# make histograms of gene expression values per sample
png(paste0("histograms_genes_per_sample_", run_id, "_part-I.png"))
hist(assay(dds_2)[,1], xlab = colnames(assay(dds_2))[1], main = "Distribution of gene expression per sample")
dev.off()

png(paste0("histograms_genes_per_sample_", run_id, "_part-II.png"))
par(mfrow=c(4,5))
hist(assay(dds_2)[,2], xlab = colnames(assay(dds_2))[2], main = NULL)
hist(assay(dds_2)[,3], xlab = colnames(assay(dds_2))[3], main = NULL)
hist(assay(dds_2)[,4], xlab = colnames(assay(dds_2))[4], main = NULL)
hist(assay(dds_2)[,5], xlab = colnames(assay(dds_2))[5], main = NULL)
hist(assay(dds_2)[,6], xlab = colnames(assay(dds_2))[6], main = NULL)
hist(assay(dds_2)[,7], xlab = colnames(assay(dds_2))[7], main = NULL)
hist(assay(dds_2)[,8], xlab = colnames(assay(dds_2))[8], main = NULL)
hist(assay(dds_2)[,9], xlab = colnames(assay(dds_2))[9], main = NULL)
hist(assay(dds_2)[,10], xlab = colnames(assay(dds_2))[10], main = NULL)
hist(assay(dds_2)[,11], xlab = colnames(assay(dds_2))[11], main = NULL)
hist(assay(dds_2)[,12], xlab = colnames(assay(dds_2))[12], main = NULL)
hist(assay(dds_2)[,13], xlab = colnames(assay(dds_2))[13], main = NULL)
hist(assay(dds_2)[,14], xlab = colnames(assay(dds_2))[14], main = NULL)
hist(assay(dds_2)[,15], xlab = colnames(assay(dds_2))[15], main = NULL)
hist(assay(dds_2)[,16], xlab = colnames(assay(dds_2))[16], main = NULL)
hist(assay(dds_2)[,17], xlab = colnames(assay(dds_2))[17], main = NULL)
hist(assay(dds_2)[,18], xlab = colnames(assay(dds_2))[18], main = NULL)
hist(assay(dds_2)[,19], xlab = colnames(assay(dds_2))[19], main = NULL)
hist(assay(dds_2)[,20], xlab = colnames(assay(dds_2))[20], main = NULL)
hist(assay(dds_2)[,21], xlab = colnames(assay(dds_2))[21], main = NULL)
dev.off()

png(paste0("histograms_genes_per_sample_", run_id, "_part-III.png"))
par(mfrow=c(4,5))
hist(assay(dds_2)[,22], xlab = colnames(assay(dds_2))[22], main = NULL)
hist(assay(dds_2)[,23], xlab = colnames(assay(dds_2))[23], main = NULL)
hist(assay(dds_2)[,24], xlab = colnames(assay(dds_2))[24], main = NULL)
hist(assay(dds_2)[,25], xlab = colnames(assay(dds_2))[25], main = NULL)
hist(assay(dds_2)[,26], xlab = colnames(assay(dds_2))[26], main = NULL)
hist(assay(dds_2)[,27], xlab = colnames(assay(dds_2))[27], main = NULL)
hist(assay(dds_2)[,28], xlab = colnames(assay(dds_2))[28], main = NULL)
hist(assay(dds_2)[,29], xlab = colnames(assay(dds_2))[29], main = NULL)
hist(assay(dds_2)[,30], xlab = colnames(assay(dds_2))[30], main = NULL)
hist(assay(dds_2)[,31], xlab = colnames(assay(dds_2))[31], main = NULL)
hist(assay(dds_2)[,32], xlab = colnames(assay(dds_2))[32], main = NULL)
hist(assay(dds_2)[,33], xlab = colnames(assay(dds_2))[33], main = NULL)
hist(assay(dds_2)[,34], xlab = colnames(assay(dds_2))[34], main = NULL)
hist(assay(dds_2)[,35], xlab = colnames(assay(dds_2))[35], main = NULL)
hist(assay(dds_2)[,36], xlab = colnames(assay(dds_2))[36], main = NULL)
hist(assay(dds_2)[,37], xlab = colnames(assay(dds_2))[37], main = NULL)
hist(assay(dds_2)[,38], xlab = colnames(assay(dds_2))[38], main = NULL)
hist(assay(dds_2)[,39], xlab = colnames(assay(dds_2))[39], main = NULL)
hist(assay(dds_2)[,40], xlab = colnames(assay(dds_2))[40], main = NULL)
hist(assay(dds_2)[,41], xlab = colnames(assay(dds_2))[41], main = NULL)
dev.off()

#############################################################################
# check the results of differential expression analysis 
# (i.e. statistical testing using the Mann-Whitney for histological and clinical features)
#############################################################################

# load the results
DE_results_dir <- paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/all_chr/mann-whitney")

# there are 3 table withg results for each feature
# => results_raw_Giant_cells_list.csv
# => results_rlog_Giant_cells_list.csv
# => results_vst_Giant_cells_list.csv

# all the features used in the analysis
# GCA_present/
# Giant_cells/
# Infiltrate_around_vasa_vasorum/
# Intima_pattern/
# Media_destruction/
# Media_pattern/
# Neoangiogenesis/
# Occlusion_grade/

# gender/
# ischaemic_features/
# jaw_claudication/
# visual_loss/

# check in Giant_cells
tab_raw_Giant_cells <- read.csv(paste0(DE_results_dir, "/Giant_cells/results_raw_Giant_cells_list.csv"))
tab_vst_Giant_cells <- read.csv(paste0(DE_results_dir, "/Giant_cells/results_vst_Giant_cells_list.csv"))
tab_rlog_Giant_cells <- read.csv(paste0(DE_results_dir, "/Giant_cells/results_rlog_Giant_cells_list.csv"))

which(tab_raw_Giant_cells$ID == "SSTR1")
which(tab_raw_Giant_cells$ID == "SSTR2")
which(tab_raw_Giant_cells$ID == "SSTR3")
which(tab_raw_Giant_cells$ID == "SSTR4")
which(tab_raw_Giant_cells$ID == "SSTR5")

which(tab_vst_Giant_cells$ID == "SSTR1")
which(tab_vst_Giant_cells$ID == "SSTR2")
which(tab_vst_Giant_cells$ID == "SSTR3")
which(tab_vst_Giant_cells$ID == "SSTR4")
which(tab_vst_Giant_cells$ID == "SSTR5")

which(tab_rlog_Giant_cells$ID == "SSTR1")
which(tab_rlog_Giant_cells$ID == "SSTR2")
which(tab_rlog_Giant_cells$ID == "SSTR3")
which(tab_rlog_Giant_cells$ID == "SSTR4")
which(tab_rlog_Giant_cells$ID == "SSTR5")


