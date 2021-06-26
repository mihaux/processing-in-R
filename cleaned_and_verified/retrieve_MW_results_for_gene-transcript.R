# this script allows to retrieve p-value and MTC_p-value for a gene or transcript of interest

# load files with results: DGE_table_<feature_name>.csv & DTE_table_<feature_name>.csv

main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES_FINAL/2021-06-06_statistical_testing_Mann-Whitney"
dir_out <- "/Users/ummz/Documents/OneDrive - University of Leeds/colaborations/colaboration_with_Jason/data_and_results"
  
setwd(dir_out)

list_all <- list.files(main_dir, pattern = "_all_results_df_", recursive = T, full.names = T)

# extract feature names to name variables
names_all_objects <- lapply(strsplit(list_all, "/"), function(x) tail(x,2)[1])

# split in 2 and add "_DGE" and "_DGE" at the end of each name
names_all_objects[seq(1,length(names_all_objects),2)] <- paste0(names_all_objects[seq(1,length(names_all_objects),2)], "_DGE")
names_all_objects[seq(2,length(names_all_objects),2)] <- paste0(names_all_objects[seq(2,length(names_all_objects),2)], "_DTE")

list_all_DGE <- list_all[which(endsWith(unlist(names_all_objects), "_DGE"))]
list_all_DTE <- list_all[which(endsWith(unlist(names_all_objects), "_DTE"))]

# load all files (it might take a while)
# DGE_all_results_df_<feature_name>.csv [use read.csv2, for ";" separated]
# DTE_all_results_df_<feature_name>.csv [use read.csv, for "," separated]

res_all_DGE <- lapply(list_all_DGE, function(x) read.csv2(x))
res_all_DTE <- lapply(list_all_DTE, function(x) read.csv(x))

# NOTE: if a gene or transcripts is not present in the results, it means that 
# it was excluded before statistical testing due to insufficient number of counts (probably just zeros everywhere)
# which(grepl("SSTR", res_all_DGE[[1]]$gene_name, fixed = TRUE)) 
# which(grepl("SSTR", res_all_DTE[[1]]$gene_name, fixed = TRUE)) 

### gene-level ###
g_sstr1 <- unique(unlist(lapply(res_all_DGE, function(x) which(x$gene_name == "SSTR1"))))
g_sstr2 <- unique(unlist(lapply(res_all_DGE, function(x) which(x$gene_name == "SSTR2"))))
#unique(unlist(lapply(res_all_DGE, function(x) which(x$gene_name == "SSTR3")))) # not present in results
#unique(unlist(lapply(res_all_DGE, function(x) which(x$gene_name == "SSTR4")))) # not present in results
g_sstr5 <- unique(unlist(lapply(res_all_DGE, function(x) which(x$gene_name == "SSTR5"))))

# check on ENSEMBL ID as well
# SSTR1   |  ENSG00000139874
#unique(unlist(lapply(res_all_DGE, function(x) which(x$gene_id == "ENSG00000139874"))))

# SSTR2   |  ENSG00000180616
#unique(unlist(lapply(res_all_DGE, function(x) which(x$gene_id == "ENSG00000180616"))))

# SSTR3   |  ENSG00000278195
#unique(unlist(lapply(res_all_DGE, function(x) which(x$gene_id == "ENSG00000278195"))))

# SSTR4   |  ENSG00000132671
#unique(unlist(lapply(res_all_DGE, function(x) which(x$gene_id == "ENSG00000132671"))))

# SSTR5   |  ENSG00000162009
#unique(unlist(lapply(res_all_DGE, function(x) which(x$gene_id == "ENSG00000162009"))))

### transcript-level ###
# "ENST00000267377.2" | SSTR1-201
tr_sstr1 <- unique(unlist(lapply(res_all_DTE, function(x) which(x$tx_id == "ENST00000267377"))))

# "ENST00000357585.3" | SSTR2-201
tr_sstr2.1 <- unique(unlist(lapply(res_all_DTE, function(x) which(x$tx_id == "ENST00000357585"))))

# "ENST00000579323.5" | SSTR2-202 # not present in results
#unique(unlist(lapply(res_all_DTE, function(x) which(x$tx_id == "ENST00000579323"))))

# "ENST00000610913.1" | SSTR3-201 # not present in results
#unique(unlist(lapply(res_all_DTE, function(x) which(x$tx_id == "ENST00000610913"))))

# "ENST00000617123.1" | SSTR3-202 # not present in results
#unique(unlist(lapply(res_all_DTE, function(x) which(x$tx_id == "ENST00000617123"))))

# "ENST00000255008.4"  | SSTR4-201 # not present in results
#unique(unlist(lapply(res_all_DTE, function(x) which(x$tx_id == "ENST00000255008"))))

# "ENST00000293897.5"  | SSTR5-201
tr_sstr5 <- unique(unlist(lapply(res_all_DTE, function(x) which(x$tx_id == "ENST00000293897"))))

# extract the results
# g_sstr1 | g_sstr2 | g_sstr5
l_g_sstr1 <- lapply(res_all_DGE, function(x) x[g_sstr1,c("log2FC", "pvalue", "qvalue")])
l_g_sstr2 <- lapply(res_all_DGE, function(x) x[g_sstr2,c("log2FC", "pvalue", "qvalue")])
l_g_sstr5 <- lapply(res_all_DGE, function(x) x[g_sstr5,c("log2FC", "pvalue", "qvalue")])

df_g_sstr1 <- data.frame(matrix(unlist(l_g_sstr1), nrow=length(l_g_sstr1), byrow=TRUE), stringsAsFactors=FALSE)
df_g_sstr2 <- data.frame(matrix(unlist(l_g_sstr2), nrow=length(l_g_sstr2), byrow=TRUE), stringsAsFactors=FALSE)
df_g_sstr5 <- data.frame(matrix(unlist(l_g_sstr5), nrow=length(l_g_sstr5), byrow=TRUE), stringsAsFactors=FALSE)

colnames(df_g_sstr1) <- c("log2FC", "pvalue", "FDR")
colnames(df_g_sstr2) <- c("log2FC", "pvalue", "FDR")
colnames(df_g_sstr5) <- c("log2FC", "pvalue", "FDR")

rownames(df_g_sstr1) <- names_all_objects[seq(1,length(names_all_objects),2)]
rownames(df_g_sstr2) <- names_all_objects[seq(1,length(names_all_objects),2)]
rownames(df_g_sstr5) <- names_all_objects[seq(1,length(names_all_objects),2)]

# save results
write.csv(df_g_sstr1, "results_Mann-Whitney_gene_SSTR1.csv")
write.csv(df_g_sstr2, "results_Mann-Whitney_gene_SSTR2.csv")
write.csv(df_g_sstr5, "results_Mann-Whitney_gene_SSTR5.csv")

# tr_sstr1 | tr_sstr2.1 | tr_sstr5
l_tr_sstr1 <- lapply(res_all_DTE, function(x) x[tr_sstr1,c("log2FC", "pvalue", "qvalue")])
l_tr_sstr2.1 <- lapply(res_all_DTE, function(x) x[tr_sstr2.1,c("log2FC", "pvalue", "qvalue")])
l_tr_sstr5 <- lapply(res_all_DTE, function(x) x[tr_sstr5,c("log2FC", "pvalue", "qvalue")])

df_tr_sstr1 <- data.frame(matrix(unlist(l_tr_sstr1), nrow=length(l_tr_sstr1), byrow=TRUE), stringsAsFactors=FALSE)
df_tr_sstr2.1 <- data.frame(matrix(unlist(l_tr_sstr2.1), nrow=length(l_tr_sstr2.1), byrow=TRUE), stringsAsFactors=FALSE)
df_tr_sstr5 <- data.frame(matrix(unlist(l_tr_sstr5), nrow=length(l_tr_sstr5), byrow=TRUE), stringsAsFactors=FALSE)

colnames(df_tr_sstr1) <- c("log2FC", "pvalue", "FDR")
colnames(df_tr_sstr2.1) <- c("log2FC", "pvalue", "FDR")
colnames(df_tr_sstr5) <- c("log2FC", "pvalue", "FDR")

rownames(df_tr_sstr1) <- names_all_objects[seq(2,length(names_all_objects),2)]
rownames(df_tr_sstr2.1) <- names_all_objects[seq(2,length(names_all_objects),2)]
rownames(df_tr_sstr5) <- names_all_objects[seq(2,length(names_all_objects),2)]

# save results
write.csv(df_tr_sstr1, "results_Mann-Whitney_transcript_SSTR1.csv")
write.csv(df_tr_sstr2.1, "results_Mann-Whitney_transcript_SSTR2.1.csv")
write.csv(df_tr_sstr5, "results_Mann-Whitney_transcript_SSTR5.csv")

################################################################################
# select only features of interests and split in clinical and histological
################################################################################

# histological features of interest
g_hist <- c("GCA_present_DGE", "Giant_cells_DGE", "Infiltrate_around_vasa_vasorum_DGE",
            "Media_destruction_DGE", "Neoangiogenesis_DGE", "Adventitia_pattern_DGE",
            "Media_pattern_DGE", "Intima_pattern_DGE", "Occlusion_grade_DGE")

tr_hist <- c("GCA_present_DTE", "Giant_cells_DTE", "Infiltrate_around_vasa_vasorum_DTE",
            "Media_destruction_DTE", "Neoangiogenesis_DTE", "Adventitia_pattern_DTE",
            "Media_pattern_DTE", "Intima_pattern_DTE", "Occlusion_grade_DTE")

# clinical features of interest
g_clinic <- c("visual_loss_DGE", "ischaemic_features_DGE", "jaw_claudication_DGE", 
              "gender_DGE", "age_DGE", "steroids_duration_DGE")

tr_clinic <- c("visual_loss_DTE", "ischaemic_features_DTE", "jaw_claudication_DTE", 
               "gender_DTE", "age_DTE", "steroids_duration_DTE")
    
### histological ###
df_g_sstr1_hist <- df_g_sstr1[which(rownames(df_g_sstr1) %in% g_hist),]
df_tr_sstr1_hist <- df_tr_sstr1[which(rownames(df_tr_sstr1) %in% tr_hist),]

df_g_sstr2_hist <- df_g_sstr2[which(rownames(df_g_sstr2) %in% g_hist),]
df_tr_sstr2.1_hist <- df_tr_sstr2.1[which(rownames(df_tr_sstr2.1) %in% tr_hist),]

df_g_sstr5_hist <- df_g_sstr5[which(rownames(df_g_sstr5) %in% g_hist),]
df_tr_sstr5_hist <- df_tr_sstr5[which(rownames(df_tr_sstr5) %in% tr_hist),]

write.csv(round(df_g_sstr1_hist, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_gene_SSTR1_selected_histo.csv")
write.csv(round(df_g_sstr2_hist, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_gene_SSTR2_selected_histo.csv")
write.csv(round(df_g_sstr5_hist, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_gene_SSTR5_selected_histo.csv")

write.csv(round(df_tr_sstr1_hist, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_transcript_SSTR1_selected_histo.csv")
write.csv(round(df_tr_sstr2.1_hist, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_transcript_SSTR2_selected_histo.csv")
write.csv(round(df_tr_sstr5_hist, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_transcript_SSTR5_selected_histo.csv")

### clinical ###
df_g_sstr1_clinic <- df_g_sstr1[which(rownames(df_g_sstr1) %in% g_clinic),]
df_tr_sstr1_clinic <- df_tr_sstr1[which(rownames(df_tr_sstr1) %in% tr_clinic),]

df_g_sstr2_clinic <- df_g_sstr2[which(rownames(df_g_sstr2) %in% g_clinic),]
df_tr_sstr2.1_clinic <- df_tr_sstr2.1[which(rownames(df_tr_sstr2.1) %in% tr_clinic),]

df_g_sstr5_clinic <- df_g_sstr5[which(rownames(df_g_sstr5) %in% g_clinic),]
df_tr_sstr5_clinic <- df_tr_sstr5[which(rownames(df_tr_sstr5) %in% tr_clinic),]

write.csv(round(df_g_sstr1_clinic, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_gene_SSTR1_selected_clinic.csv")
write.csv(round(df_g_sstr2_clinic, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_gene_SSTR2_selected_clinic.csv")
write.csv(round(df_g_sstr5_clinic, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_gene_SSTR5_selected_clinic.csv")

write.csv(round(df_tr_sstr1_clinic, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_transcript_SSTR1_selected_clinic.csv")
write.csv(round(df_tr_sstr2.1_clinic, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_transcript_SSTR2_selected_clinic.csv")
write.csv(round(df_tr_sstr5_clinic, digits = 5), "results_Mann-Whitney_selected/results_Mann-Whitney_transcript_SSTR5_selected_clinic.csv")
















