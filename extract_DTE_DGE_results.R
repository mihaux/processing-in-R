# script to extract all the DTE and DGE results and create summary table

main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES_archived/swish_v2/NO_OUTLIERS_40"

list_all_nb_down_up <- list.files(file.path(main_dir, "output_all_NO-OUTLIERS"), pattern = "nb_down-up", recursive = T, full.names = T)
list_all_nb_significant <- list.files(file.path(main_dir, "output_all_NO-OUTLIERS"), pattern = "nb_significant", recursive = T, full.names = T)

list_1_22_nb_down_up <- list.files(file.path(main_dir, "output_1-22_NO-OUTLIERS"), pattern = "nb_down-up", recursive = T, full.names = T)
list_1_22_nb_significant <- list.files(file.path(main_dir, "output_1-22_NO-OUTLIERS"), pattern = "nb_significant", recursive = T, full.names = T)

list_xy_nb_down_up <- list.files(file.path(main_dir, "output_X-Y_NO-OUTLIERS"), pattern = "nb_down-up", recursive = T, full.names = T)
list_xy_nb_significant <- list.files(file.path(main_dir, "output_X-Y_NO-OUTLIERS"), pattern = "nb_significant", recursive = T, full.names = T)

# extract feture names to name variables
names_all_objects <- sapply(strsplit(list.files(file.path(main_dir, "output_all_NO-OUTLIERS"), pattern = "nb_down-up", recursive = T) , "/"), "[", 1)
# SAME FOR: names_all_nb_significant, names_1_22_nb_down_up, names_1_22_nb_significant, names_xy_nb_down_up, names_xy_nb_significant
names_all_objects[seq(1,length(names_all_objects),2)] <- paste0(names_all_objects[seq(1,length(names_all_objects),2)], "_DGE")
names_all_objects[seq(2,length(names_all_objects),2)] <- paste0(names_all_objects[seq(2,length(names_all_objects),2)], "_DTE")

# there are 4 files to be processed foe each feature (2 for DTE & 2 for DGE):
# examples:
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
res_1_22_nb_down_up <- lapply(list_1_22_nb_down_up, function(x) scan(x, what = "character"))
res_xy_nb_down_up <- lapply(list_xy_nb_down_up, function(x) scan(x, what = "character"))

names(res_all_nb_down_up) <- names_all_objects
names(res_1_22_nb_down_up) <- names_all_objects
names(res_xy_nb_down_up) <- names_all_objects

# check if there are statistically signifant results for each feature
str(res_all_nb_down_up)   # Infiltrate_around_vasa_vasorum_DTE, Neoangiogenesis_DGE, Neoangiogenesis_DTE
str(res_1_22_nb_down_up)  # Infiltrate_around_vasa_vasorum_DTE, Neoangiogenesis_DGE, Neoangiogenesis_DTE
str(res_xy_nb_down_up)    # Infiltrate_around_vasa_vasorum_DTE, Neoangiogenesis_DGE, Neoangiogenesis_DTE

# extract lines 26-29 from each
res_all_nb_down_up_extracted <- lapply(res_all_nb_down_up, function(x) x[c(27,29)])
res_1_22_nb_down_up_extracted <- lapply(res_1_22_nb_down_up, function(x) x[c(27,29)])
res_xy_nb_down_up_extracted <- lapply(res_xy_nb_down_up, function(x) x[c(27,29)])

# save tables with results
setwd("/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES_archived/swish_v2/NO_OUTLIERS_40/results_extracted")

# ----- all -----
tb_all_nb_down_up <- data.frame(negative_FC=unlist(lapply(res_all_nb_down_up_extracted, function(x) x[1])),
                                positive_FC=unlist(lapply(res_all_nb_down_up_extracted, function(x) x[2])))

tb_all_nb_down_up_DGE <- tb_all_nb_down_up[seq(1,nrow(tb_all_nb_down_up),2),]
tb_all_nb_down_up_DTE <- tb_all_nb_down_up[seq(2,nrow(tb_all_nb_down_up),2),]

write.csv(tb_all_nb_down_up_DGE, "table_all_nb_down_up_DGE.csv")
write.csv(tb_all_nb_down_up_DTE, "table_all_nb_down_up_DTE.csv")

# ----- 1-22 -----
tb_1_22_nb_down_up <- data.frame(negative_FC=unlist(lapply(res_1_22_nb_down_up_extracted, function(x) x[1])),
                                positive_FC=unlist(lapply(res_1_22_nb_down_up_extracted, function(x) x[2])))

tb_1_22_nb_down_up_DGE <- tb_1_22_nb_down_up[seq(1,nrow(tb_1_22_nb_down_up),2),]
tb_1_22_nb_down_up_DTE <- tb_1_22_nb_down_up[seq(2,nrow(tb_1_22_nb_down_up),2),]

write.csv(tb_1_22_nb_down_up_DGE, "table_1_22_nb_down_up_DGE.csv")
write.csv(tb_1_22_nb_down_up_DTE, "table_1_22_nb_down_up_DTE.csv")

# ----- X-Y -----
tb_xy_nb_down_up <- data.frame(negative_FC=unlist(lapply(res_xy_nb_down_up_extracted, function(x) x[1])),
                                 positive_FC=unlist(lapply(res_xy_nb_down_up_extracted, function(x) x[2])))

tb_xy_nb_down_up_DGE <- tb_xy_nb_down_up[seq(1,nrow(tb_xy_nb_down_up),2),]
tb_xy_nb_down_up_DTE <- tb_xy_nb_down_up[seq(2,nrow(tb_xy_nb_down_up),2),]

write.csv(tb_xy_nb_down_up_DGE, "table_xy_nb_down_up_DGE.csv")
write.csv(tb_xy_nb_down_up_DTE, "table_xy_nb_down_up_DTE.csv")

######## in DXE_nb_significant_XXX.txt: ######## 
# line 10: number of transcript/genes with p-values lower than 5%
# line 19: number of transcript/genes with p-values lower than 1%

res_all_nb_significant <- lapply(list_all_nb_significant, function(x) scan(x, what = "character"))
res_1_22_nb_significant <- lapply(list_1_22_nb_significant, function(x) scan(x, what = "character"))
res_xy_nb_significant <- lapply(list_xy_nb_significant, function(x) scan(x, what = "character"))

names(res_all_nb_significant) <- names_all_objects
names(res_1_22_nb_significant) <- names_all_objects
names(res_xy_nb_significant) <- names_all_objects
  
lapply(res_all_nb_significant, function(x) x[8:length(x)])
lapply(res_1_22_nb_significant, function(x) x[8:length(x)])
lapply(res_xy_nb_significant, function(x) x[8:length(x)])

# ----- all -----

# ----- 1-22 -----

# ----- X-Y -----
