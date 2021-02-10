# make plots of expression profile of a selected transcript (using results from Salmon)
# data types: Raw 

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
data_dir <- paste0(main_dir,"/ANALYSES_archived/run_13_Jan21/processed/")

# define directory for results (OUTPUT)
dir_out <- paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/expression_profiles")
setwd(dir_out)

# load data 
df <- read.csv(paste0(data_dir, "counts_transcript_level.csv"), row.names = 1)

# define running ID (either "raw", "vst" pr "rlog")
run_id <- "raw"

#############################################################################
# => in df, look for:
# SSTR1-201   |  ENST00000267377.3
df_1 <- which(startsWith(rownames(df), "ENST00000267377"))    # 86319
rownames(df)[df_1] # "ENST00000267377.2"

# SSTR2-201   |  ENST00000357585.4
# SSTR2-202   |  ENST00000579323.5
df_2.1 <- which(startsWith(rownames(df), "ENST00000357585"))  # 135293
df_2.2 <- which(startsWith(rownames(df), "ENST00000579323"))  # 135292
rownames(df)[df_2.1]  # "ENST00000357585.3"
rownames(df)[df_2.2]  # "ENST00000579323.5"

# SSTR3-201   |  ENST00000610913.2
# SSTR3-202   |  ENST00000617123.1
df_3.1 <- which(startsWith(rownames(df), "ENST00000610913"))  # 127236
df_3.2 <- which(startsWith(rownames(df), "ENST00000617123"))  # 127237
rownames(df)[df_3.1]  # "ENST00000610913.1"
rownames(df)[df_3.2]  # "ENST00000617123.1"

# SSTR4-201   |  ENST00000255008.5
df_4 <- which(startsWith(rownames(df), "ENST00000255008"))    # 144722
rownames(df)[df_4] # "ENST00000255008.4"

# SSTR5-201   |  ENST00000293897.5
df_5 <- which(startsWith(rownames(df), "ENST00000293897"))    # 106422
rownames(df)[df_5] # "ENST00000293897.5"

# get min, max, mean and median for each sample
df_min <- apply(df, 2, min)
df_max <- apply(df, 2, max)
df_mean <- apply(df, 2, mean)
df_median <- apply(df, 2, median)

# save results table
res_tab_1 <- data.frame(min=df_min, max=df_max, mean=df_mean, median=df_median)
if(output_save==TRUE){ write.csv(res_tab_1, file = paste0("table_transcripts_", run_id,".csv")) }

# create a matrix with all transcripts of interest
mat_1 <- rbind(df[df_1,], 
               df[df_2.1,], df[df_2.2,], 
               df[df_3.1,], df[df_3.2,], 
               df[df_4,],  
               df[df_5,], 
               df_max,
               df_mean)

rownames(mat_1) <- c("SSTR1-201", "SSTR2-201", "SSTR2-202", "SSTR3-201", "SSTR3-202", "SSTR4-201", "SSTR5-201", "max", "mean")
  
# write a table 
write.csv(mat_1, "transcripts_of_interest-counts.csv")



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


