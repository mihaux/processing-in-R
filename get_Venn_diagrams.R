# this script creates Venn diagrams for transcript lists (output from test_Mann-Whitney.R)

# install (if necessary) and load package
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
suppressMessages(library(VennDiagram))
suppressMessages(library(stringr))

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

#args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=2) {
  stop("2 arguments must be supplied: 
       \n(1 - input) path to the directory with transcript list data and, 
       \n(2 - output) path where output files should be stored", call.=FALSE)
}

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3

args <- c(paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/"),
          paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/venn_diagrams/"))

# for SE
# /ANALYSES/run_12_Aug20/6_downstream/SE/DESeq2_analysis/all_chr/mann-whitney
# /ANALYSES/run_12_Aug20/6_downstream/SE/DESeq2_analysis/no_chrXY/mann-whitney

# for PE
# /ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/all_chr/mann-whitney
# /ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/no_chrXY/mann-whitney

# (OUT) /ANALYSES/run_12_Aug20/6_downstream/venn_diagrams

# NOTE: there are 3 files for each feature
# => results_raw_GCA_present_list.csv
# => results_rlog_GCA_present_list.csv
# => results_vst_GCA_present_list.csv

# Example of usage: 
# Rscript 

# clinical features
# => gender
# => visual_loss
# => jaw_claudication
# => ischaemic_features

# histological features
# => GCA_present
# => Giant_cells
# => Media_destruction
# => Occlusion_grade
# => Neoangiogenesis
# => Intima_pattern
# => Infiltrate_around_vasa_vasorum

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[2], sep="\n")
setwd(args[2])

# load count data (cols -> samples | rows -> genes)
df <- read.csv(args[1], row.names = 1, header = TRUE)     # data.frame with counts only

# get all files to be loaded
input_files <- list.files(path = args[1], recursive = TRUE, pattern = "list.csv", full.names = TRUE)

# load all files
myfiles = lapply(input_files, read.csv)

# extract names for each file
extracted_names <- str_split(input_files, "/", n = Inf, simplify = FALSE)

# create final names and add it to the list object
names(myfiles) <- lapply(extracted_names, function(x) paste0(x[10], "_", x[12], "_", x[15]))

# OVERLAPS TO BE CHECKED:
# (1) for each feature => between the 3 methods
# (2) between all_chr and no_chrXY
# (3) for each feature => between SE and PE
# (4) between different sets of features (if more than 5 then just get intersections of all)



# save transformed data and unnormalised as well
#save(dds_all, file="Raw_DESeq_dataset_all.Rda")
#save(vst_all, file="Normalised_DESeq_vst_dataset_all.Rda")
#save(rlog_all, file="Normalised_DESeq_rlog_dataset_all.Rda")

#save(dds_all, file="Raw_DESeq_dataset_noXY.Rda")
#save(vst_all, file="Normalised_DESeq_vst_dataset_noXY.Rda")
#save(rlog_all, file="Normalised_DESeq_rlog_dataset_noXY.Rda")


