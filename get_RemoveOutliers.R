# date: 2021-01-04
# this script removes outliers and save count matrices as .csv files as well as spreadsheet with metadata

library(stringr)

# get working directory to recognise the machine
w_dir <- getwd()

# create a shortcut for the OneDrive directory where all files are stored
if(startsWith(w_dir, "/Users/michal")){           
  main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"    # on my mac
} else if (startsWith(w_dir, "/Users/ummz")) {    
  main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds"                # on uni mac    
} else {
  print("Unrecognised machine.")
}

outliers <- "ID-8546"

# define input directories 
dir_counts_raw <- paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/all_chr/INPUT_counts/Raw_DESeq_dataset_all.Rda")
dir_counts_vst <- paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/all_chr/INPUT_counts/Normalised_DESeq_vst_dataset_all.Rda")
dir_counts_rlog <- paste0(main_dir,"/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/all_chr/INPUT_counts/Normalised_DESeq_rlog_dataset_all.Rda")
dir_metadata <- paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv")
dir_histodata <- paste0(main_dir, "/data/metadata/slide_scores/slide_scores_v6.csv")

# load counts 
load(dir_counts_raw)    # dds_all
load(dir_counts_vst)    # vst_all
load(dir_counts_rlog)   # rlog_all

# change data format to matrix
dat_raw <- as.matrix(assay(dds_all))
dat_vst <- as.matrix(assay(vst_all))
dat_rlog <- as.matrix(assay(rlog_all))

# add "ID_" to "ID-" in all count matrices
colnames(dat_raw) <- str_replace_all(colnames(dat_raw), "ID_", "ID-")
colnames(dat_vst) <- str_replace_all(colnames(dat_vst), "ID_", "ID-")
colnames(dat_rlog) <- str_replace_all(colnames(dat_rlog), "ID_", "ID-")

# remove outlier sample(s)
dat_raw_no_outliers <- dat_raw[,-which(colnames(dat_raw) == outliers)]
dat_vst_no_outliers <- dat_vst[,-which(colnames(dat_vst) == outliers)]
dat_rlog_no_outliers <- dat_rlog[,-which(colnames(dat_rlog) == outliers)]

# save modified metadata 
write.csv(dat_raw_no_outliers, file="/Users/ummz/Documents/OneDrive - University of Leeds/data/count_matrices/outliers_excluded/counts_raw_no_outliers.csv")
write.csv(dat_vst_no_outliers, file="/Users/ummz/Documents/OneDrive - University of Leeds/data/count_matrices/outliers_excluded/counts_vst_no_outliers.csv")
write.csv(dat_rlog_no_outliers, file="/Users/ummz/Documents/OneDrive - University of Leeds/data/count_matrices/outliers_excluded/counts_rlog_no_outliers.csv")

# save unmodified data as well (TO BE REMOVED LATER)
write.csv(dat_raw, file="/Users/ummz/Documents/OneDrive - University of Leeds/data/count_matrices/full/counts_raw.csv")
write.csv(dat_vst, file="/Users/ummz/Documents/OneDrive - University of Leeds/data/count_matrices/full/counts_vst.csv")
write.csv(dat_rlog, file="/Users/ummz/Documents/OneDrive - University of Leeds/data/count_matrices/full/counts_rlog.csv")

#--------------------------------------------------------------------------------#

# load clinical data 
df_meta <- read.csv(dir_metadata, row.names = 1, header = TRUE)

# add "ID-" to all rownames
rownames(df_meta) <- paste0("ID-", rownames(df_meta))

# remove outlier sample (ID-8546)
df_meta_no_outliers <- df_meta[-which(rownames(df_meta) == outliers),]

# save modified metadata 
write.csv(df_meta_no_outliers, file="/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/outliers_excluded/cic_clinical_data_v2_summary_ORDERED_outliers_excluded.csv")


#--------------------------------------------------------------------------------#

# load histological data 
df_histo <- read.csv(dir_histodata, row.names = 1, header = TRUE)

# add "ID-" to all rownames
rownames(df_histo) <- paste0("ID-", rownames(df_histo))

# remove outlier sample (ID-8546)
df_histo_no_outliers <- df_histo[-which(rownames(df_histo) == outliers),]

# save modified metadata 
write.csv(df_histo_no_outliers, file="/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/outliers_excluded/slide_scores_v6_outliers_excluded.csv")





