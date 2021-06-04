# overlapping results from Swish analysis

library(VennDiagram)
library(RColorBrewer)

# load the following output files

main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES_archived/swish_v2/NO_OUTLIERS_40"

# NOTE: use only the results from output_all_NO-OUTLIERS/
# ! could check for others later: output_1-22_NO-OUTLIERS/ & output_X-Y_NO-OUTLIERS/

# Adventitia_pattern/
# GCA_present/
# Giant_cells/
# Infiltrate_around_vasa_vasorum/
# Intima_pattern/
# Media_destruction/
# Media_pattern/
# Neoangiogenesis/
# Occlusion_grade/

# need to extract pvalue and qvalue colummns from these files
# DGE_significant_list_Adventitia_pattern.csv 
# DTE_significant_list_Adventitia_pattern.csv

files_in <- list.files(file.path(main_dir, "output_all_NO-OUTLIERS"), recursive = T, full.names = T, pattern = "_significant_list_")

# extract feture names to name variables
names_all_objects <- sapply(strsplit(list.files(file.path(main_dir, "output_all_NO-OUTLIERS"), pattern = "_significant_list_", recursive = T) , "/"), "[", 1)

# these 2 files could also be useful
#DGE_table_negative_log_fold_change_Adventitia_pattern.csv
#DGE_table_positive_log_fold_change_Adventitia_pattern.csv

#DTE_table_negative_log_fold_change_Adventitia_pattern.csv
#DTE_table_positive_log_fold_change_Adventitia_pattern.csv
  
# load all .csv files
res_all_significant_list <- lapply(files_in, function(x) read.csv(x))
names(res_all_significant_list) <- names_all_objects

# split in 'DTE' and 'DGE'
res_all_significant_list_DGE <- res_all_significant_list[seq(1, length(res_all_significant_list), by=2)]
res_all_significant_list_DTE <- res_all_significant_list[seq(2, length(res_all_significant_list), by=2)]

# filter the results to get only qvalues > 0.01
# >>> gene-level <<<
res_all_significant_list_DGE_01_ind <- lapply(res_all_significant_list_DGE, function(x) which(x$qvalue < 0.01))
# NOTE: no results for Infiltrate_around_vasa_vasorum (4), Media_destruction (6), Neoangiogenesis (8) and Occlusion_grade (9)

final_list_DGE_01 <- list(as.character(res_all_significant_list_DGE[[1]][res_all_significant_list_DGE_01_ind[[1]],]$gene_name),
                       as.character(res_all_significant_list_DGE[[2]][res_all_significant_list_DGE_01_ind[[2]],]$gene_name),
                       as.character(res_all_significant_list_DGE[[3]][res_all_significant_list_DGE_01_ind[[3]],]$gene_name),
                       as.character(res_all_significant_list_DGE[[5]][res_all_significant_list_DGE_01_ind[[5]],]$gene_name),
                       as.character(res_all_significant_list_DGE[[7]][res_all_significant_list_DGE_01_ind[[7]],]$gene_name))

names(final_list_DGE_01) <- names(res_all_significant_list_DGE)[-c(4,6,8,9)]

# get intersection across all
Reduce(intersect, final_list_DGE_01)

# get intersection across all, but Media_pattern
Reduce(intersect, final_list_DGE_01[-5])

# get intersection across all, but  GCA_present and Media_pattern
Reduce(intersect, final_list_DGE_01[-c(2,5)])

myCol <- brewer.pal(5, "Pastel2")

# create a venn diagrams (5 elements max)
venn.diagram(x = final_list_DGE_01, category.names = names(final_list_DGE_01),
             filename = 'venn_diagramm_DGE_01_all.png',output=TRUE, fill = myCol, cex = 1.3)

venn.diagram(x = final_list_DGE_01[-5], category.names = names(final_list_DGE_01[-5]),
             filename = 'venn_diagramm_DGE_01_without_1.png',output=TRUE, fill = myCol[-5], cex = 1.3)

venn.diagram(x = final_list_DGE_01[-c(2,5)], category.names = names(final_list_DGE_01[-c(2,5)]),
             filename = 'venn_diagramm_DGE_01_without_2.png',output=TRUE, fill = myCol[-c(2,5)], cex = 1.3)

# >>> transcript-level <<<
res_all_significant_list_DTE_01_ind <- lapply(res_all_significant_list_DTE, function(x) which(x$qvalue < 0.01))
# NOTE: no results for Infiltrate_around_vasa_vasorum (4), Neoangiogenesis (8) and Occlusion_grade (9)

final_list_DTE_01 <- list(as.character(res_all_significant_list_DTE[[1]][res_all_significant_list_DTE_01_ind[[1]],]$X),
                          as.character(res_all_significant_list_DTE[[2]][res_all_significant_list_DTE_01_ind[[2]],]$X),
                          as.character(res_all_significant_list_DTE[[3]][res_all_significant_list_DTE_01_ind[[3]],]$X),
                          as.character(res_all_significant_list_DTE[[5]][res_all_significant_list_DTE_01_ind[[5]],]$X),
                          as.character(res_all_significant_list_DTE[[6]][res_all_significant_list_DTE_01_ind[[6]],]$X),
                          as.character(res_all_significant_list_DTE[[7]][res_all_significant_list_DTE_01_ind[[7]],]$X))

names(final_list_DTE_01) <- names(res_all_significant_list_DTE)[-c(4,8,9)]

# get intersection across all
Reduce(intersect, final_list_DTE_01)

#
Reduce(intersect, final_list_DTE_01[-5])

#
Reduce(intersect, final_list_DTE_01[-c(5,6)])

# get intersections of $X column (i.e. ENSEMBL IDs)

# create a venn diagrams (5 elements max)
venn.diagram(x = final_list_DTE_01[-5], category.names = names(final_list_DTE_01[-5]),
  filename = 'venn_diagramm_DTE_01_without_1.png',output=TRUE, fill = myCol, cex = 1.3)

venn.diagram(x = final_list_DTE_01[-c(5,6)], category.names = names(final_list_DTE_01[-c(5,6)]),
             filename = 'venn_diagramm_DTE_01_without_2.png',output=TRUE, fill = myCol[-c(5)], cex = 1.3)

venn.diagram(x = final_list_DTE_01[-c(2,5,6)], category.names = names(final_list_DTE_01[-c(2,5,6)]),
             filename = 'venn_diagramm_DTE_01_without_3.png',output=TRUE, fill = myCol[-c(2,5)], cex = 1.3)




  