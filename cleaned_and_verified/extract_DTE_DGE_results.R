# script to extract all the DTE and DGE results and create summary table

# UPDATED on 2021-06-04
# changes log:
# - removed repeated code, it can be run for one output directory at once 
# (i.e. when running for all, chr1-22 and chrX-Y, it needs to be launched three times separately for each folder)

# set of output files (per feature)
# NOTE: there's exactly same set of files starting with "DTE" instead of "DGE"  (18 files per directory)

# DGE_all_results_df_<feature_name>.csv
# DGE_histogram_pvalues_<feature_name>.png                => image
# DGE_MA_plot_<feature_name>.png                          => image
# DGE_nb_down-up_<feature_name>.txt                       => used in this file
# DGE_nb_significant_<feature_name>.txt                   => used in this file
# DGE_significant_list_<feature_name>.csv
# DGE_significant_list_pval_<feature_name>.csv
# DGE_significant_list_qval_<feature_name>.csv
# DGE_table_<feature_name>.csv
# DGE_table_negative_log_fold_change_<feature_name>.csv
# DGE_table_positive_log_fold_change_<feature_name>.csv

main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES_FINAL/2021-06-06_statistical_testing_Mann-Whitney"

# name of the run (usually date or some other specification)
run_ID <- "2021-06-06"

list_all_nb_down_up <- list.files(main_dir, pattern = "nb_down-up", recursive = T, full.names = T)
list_all_nb_significant <- list.files(main_dir, pattern = "nb_significant", recursive = T, full.names = T)

# extract feature names to name variables
names_all_objects <- sapply(strsplit(list.files(main_dir, pattern = "nb_down-up", recursive = T) , "/"), "[", 1)

# add "_DGE" and "_DGE" at the end of each name
names_all_objects[seq(1,length(names_all_objects),2)] <- paste0(names_all_objects[seq(1,length(names_all_objects),2)], "_DGE")
names_all_objects[seq(2,length(names_all_objects),2)] <- paste0(names_all_objects[seq(2,length(names_all_objects),2)], "_DTE")

# there are 4 files to be processed foe each feature (2 for DTE & 2 for DGE); EXAMPLES:
# DGE_nb_down-up_Adventitia_pattern.txt & DGE_nb_significant_Adventitia_pattern.txt
# DTE_nb_down-up_Adventitia_pattern.txt & DTE_nb_significant_Adventitia_pattern.txt

######## in DXE_nb_down-up_XXX.txt ######## 

# total length of the results vector: 86 689
# NOTE: if(length == 25) => no statistically significant results
# line 26: "TRUE"
# line 27: number of transcript/genes with p-values lower than 5% with negative fold change
# line 28: number of transcript/genes with p-values lower than 5% with fold change = 0
# line 29: number of transcript/genes with p-values lower than 5% with positive fold change

res_all_nb_down_up <- lapply(list_all_nb_down_up, function(x) scan(x, what = "character"))

names(res_all_nb_down_up) <- names_all_objects

# check if there are statistically signifant results for each feature
#str(res_all_nb_down_up)  

# NOTE: normally there should be 29 elements for each 
# print those features that don't have 29 elements
names(res_all_nb_down_up)[which(unlist(lapply(res_all_nb_down_up, function(x) length(x))) != 29)]

# "ischaemic_features_DGE" "ischaemic_features_DTE" "jaw_claudication_DGE" "jaw_claudication_DTE"   "steroids_duration_DGE"  "steroids_duration_DTE" 

# extract lines 26-29 from each
res_all_nb_down_up_extracted <- lapply(res_all_nb_down_up, function(x) x[c(27,29)])

# save tables with results
setwd(file.path(main_dir, "x_EXTRACTED"))

tb_all_nb_down_up <- data.frame(negative_FC=unlist(lapply(res_all_nb_down_up_extracted, function(x) x[1])),
                                positive_FC=unlist(lapply(res_all_nb_down_up_extracted, function(x) x[2])))

#------------------------------------------------------------------------------------------#
# replace "NA" with zeros (it won't be needed later)
tb_all_nb_down_up_mod <- tb_all_nb_down_up
tb_all_nb_down_up_mod[which(is.na(tb_all_nb_down_up_mod[,1])),1] <- 0
tb_all_nb_down_up_mod[which(is.na(tb_all_nb_down_up_mod[,2])),2] <- 0

tb_all_nb_down_up_DGE_mod <- tb_all_nb_down_up_mod[seq(1,nrow(tb_all_nb_down_up_mod),2),]
tb_all_nb_down_up_DTE_mod <- tb_all_nb_down_up_mod[seq(2,nrow(tb_all_nb_down_up_mod),2),]

write.csv(tb_all_nb_down_up_DGE_mod, paste0(run_ID, "_DGE_summary_nb_down_up_mod.csv"))
write.csv(tb_all_nb_down_up_DTE_mod, paste0(run_ID, "_DTE_summary_nb_down_up_mod.csv"))
#------------------------------------------------------------------------------------------#

tb_all_nb_down_up_DGE <- tb_all_nb_down_up[seq(1,nrow(tb_all_nb_down_up),2),]
tb_all_nb_down_up_DTE <- tb_all_nb_down_up[seq(2,nrow(tb_all_nb_down_up),2),]

write.csv(tb_all_nb_down_up_DGE, paste0(run_ID, "_DGE_summary_nb_down_up.csv"))
write.csv(tb_all_nb_down_up_DTE, paste0(run_ID, "_DTE_summary_nb_down_up.csv"))

######## in DXE_nb_significant_XXX.txt: ######## 
# line 10: number of transcript/genes with p-values lower than 5%
# line 19: number of transcript/genes with p-values lower than 1%

res_all_nb_significant <- lapply(list_all_nb_significant, function(x) scan(x, what = "character"))

names(res_all_nb_significant) <- names_all_objects

# split in 2 to separate DGE and DTE
res_DGE_nb_significant <- res_all_nb_significant[seq(1,length(res_all_nb_significant),2)]
res_DTE_nb_significant <- res_all_nb_significant[seq(2,length(res_all_nb_significant),2)]

# 6 | 18        - p-value
# 7 | 9         - <
# 8 - 0.05 | 20 - 0.01
# 9 | 21        - "FALSE"
# 10 | 22       - "TRUE"
# 11 | 23       - nb of false
# 12 | 24       - nb of true

# 30 | 42       - q-value
# 31 | 43       - <
# 32 - 0.05 | 44 - 0.01
# 33 | 45       - "FALSE"
# 34 | 46       - "TRUE"
# 35 | 47       - nb of false
# 36 | 48       - nb of true

# gene-level
df_res_DGE <- data.frame(features=names(res_DGE_nb_significant),
                         p_val.05=unlist(lapply(res_DGE_nb_significant, function(x) x[12])),
                         p_val.01=unlist(lapply(res_DGE_nb_significant, function(x) x[24])),
                         q_val.05=unlist(lapply(res_DGE_nb_significant, function(x) x[36])),
                         q_val.01=unlist(lapply(res_DGE_nb_significant, function(x) x[48])),
                         row.names = NULL)

write.csv(df_res_DGE, paste0(run_ID, "_DGE_summary_nb_significant.csv"))

# transcript-level
df_res_DTE <- data.frame(features=names(res_DTE_nb_significant),
                         p_val.05=unlist(lapply(res_DTE_nb_significant, function(x) x[12])),
                         p_val.01=unlist(lapply(res_DTE_nb_significant, function(x) x[24])),
                         q_val.05=unlist(lapply(res_DTE_nb_significant, function(x) x[36])),
                         q_val.01=unlist(lapply(res_DTE_nb_significant, function(x) x[48])),
                         row.names = NULL)

write.csv(df_res_DTE, paste0(run_ID, "_DTE_summary_nb_significant.csv"))

# get total number of genes / transcripts
sink(file = paste0(run_ID, "_total_nb_of_genes_and_transcripts.txt"))

cat("Total number of genes: ")
cat(unique(as.numeric(unlist(lapply(res_DGE_nb_significant, function(x) x[11]))) + as.numeric(unlist(lapply(res_DGE_nb_significant, function(x) x[12])))))
cat("\n")
cat("Total number of transcripts: ")
cat(unique(as.numeric(unlist(lapply(res_DTE_nb_significant, function(x) x[11]))) + as.numeric(unlist(lapply(res_DTE_nb_significant, function(x) x[12])))))

sink()





