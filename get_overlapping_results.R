# load all the results from statistical testing (on gene-level) and get overlapping list

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
data_dir <- paste0(main_dir,"/ANALYSES_archived/run_13_Jan21/statistical_testing/output/")

# define directory for results (OUTPUT)
dir_out <- paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/statistical_testing/output/overlapping/")
setwd(dir_out)

# define paths to the results
res_dir_clinical <- list.files(paste0(data_dir, "clinical"), full.names = TRUE)
res_dir_histological <- list.files(paste0(data_dir, "histological"), full.names = TRUE)

# load result files sorted by p-adj (raw and VST for each feature)
# table_sorted_by_padjusted_<feature>_raw_gene-level.csv
# table_sorted_by_padjusted_<feature>_VST_gene-level.csv

# get list of all results files
histo_list <- list.files(res_dir_histological, pattern = "padjusted", full.names = TRUE)

# filter to get only "transcript-level.csv" (keep every second)
histo_list_fin <- histo_list[seq(1, length(histo_list), by = 2)]

histo_list_files <- lapply(histo_list_fin, function(x) read.csv2(x, sep = ";", header = TRUE, stringsAsFactors = FALSE))

names_1 <- list.files(res_dir_histological, pattern = "padjusted")
names_2 <- names_1[seq(1, length(histo_list), by = 2)]
names_3 <- substring(names_2, first = 27)
names_4 <- stringi::stri_reverse(names_3)
names_5 <- substring(names_4, first = 16)
names_6 <- stringi::stri_reverse(names_5)

names(histo_list_files) <- names_6

# separate raw and VST results
histo_list_files_raw <- histo_list_files[seq(1, length(histo_list_files), by=2)]
histo_list_files_VST <- histo_list_files[seq(2, length(histo_list_files), by=2)]
  
# get the number of statisticaly significant results
table_res_fdr <- data.frame(feature=names(histo_list_files_raw),
                            raw_data=c(length(which(histo_list_files_raw[[1]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[2]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[3]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[4]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[5]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[6]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[7]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[8]]$fdr.pvalue < 0.05))),
                            vst_data=c(length(which(histo_list_files_VST[[1]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[2]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[3]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[4]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[5]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[6]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[7]]$fdr.pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[8]]$fdr.pvalue < 0.05))))

table_res_pval <- data.frame(feature=names(histo_list_files_raw),
                            raw_data=c(length(which(histo_list_files_raw[[1]]$pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[2]]$pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[3]]$pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[4]]$pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[5]]$pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[6]]$pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[7]]$pvalue < 0.05)),
                                       length(which(histo_list_files_raw[[8]]$pvalue < 0.05))),
                            vst_data=c(length(which(histo_list_files_VST[[1]]$pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[2]]$pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[3]]$pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[4]]$pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[5]]$pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[6]]$pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[7]]$pvalue < 0.05)),
                                       length(which(histo_list_files_VST[[8]]$pvalue < 0.05))))




