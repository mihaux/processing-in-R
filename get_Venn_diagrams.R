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
#df <- read.csv(args[1], row.names = 1, header = TRUE)     # data.frame with counts only

# get all files to be loaded
input_files <- list.files(path = args[1], recursive = TRUE, pattern = "list.csv", full.names = TRUE)

# load all files
myfiles = lapply(input_files, read.csv)

# extract names for each file
extracted_names <- str_split(input_files, "/", n = Inf, simplify = FALSE)

# create final names and add it to the list object
names(myfiles) <- lapply(extracted_names, function(x) paste0(x[10], "_", x[12], "_", x[15]))

# OVERLAPS TO BE CHECKED:
# (1) for each feature => between the 2 normalisation methods
# (2) between all_chr and no_chrXY
# (3) for each feature => between SE and PE
# (4) between different sets of features (if more than 5 then just get intersections of all)

# PE => 1-66
#   => all_chr  -> 1 - 33
#   => no_chrXY -> 34 - 66
myfiles_PE <- myfiles[1:66] ; extracted_names_PE <- extracted_names[1:66]
myfiles_PE_all_chr <- myfiles_PE[1:33]; extracted_names_PE_all_chr <- extracted_names_PE[1:33]
myfiles_PE_no_chrXY <- myfiles_PE[34:66]; extracted_names_PE_no_chrXY <- extracted_names_PE[34:66]

# SE => 67-132
#   => all_chr  -> 1 - 33 (67 - 99)
#   => no_chrXY -> 34 - 66 (100 - 132)
myfiles_SE <- myfiles[67:132] ; extracted_names_SE <- extracted_names[67:132]
myfiles_SE_all_chr <- myfiles_SE[1:33]; extracted_names_SE_all_chr <- extracted_names_SE[1:33]
myfiles_SE_no_chrXY <- myfiles_SE[34:66]; extracted_names_SE_no_chrXY <- extracted_names_SE[34:66]

#substr(names(myfiles_SE[1:33]), 20, 63)
#substr(names(myfiles_SE[34:66]), 20, 63)
# NOTE: in each trio: raw | rlog | vst => so 2nd and 3rd to be compared

# => gender                           => 4:6
int_SE_gender_all_chr <- intersect(as.character(myfiles_SE_all_chr[[5]]$ID), as.character(myfiles_SE_all_chr[[6]]$ID))
int_SE_gender_no_chrXY <- intersect(as.character(myfiles_SE_no_chrXY[[5]]$ID), as.character(myfiles_SE_no_chrXY[[6]]$ID))
int_PE_gender_all_chr <- intersect(as.character(myfiles_PE_all_chr[[5]]$ID), as.character(myfiles_PE_all_chr[[6]]$ID))
int_PE_gender_no_chrXY <- intersect(as.character(myfiles_PE_no_chrXY[[5]]$ID), as.character(myfiles_PE_no_chrXY[[6]]$ID))

# => visual_loss                      => 31:33
int_SE_visual_loss_all_chr <- intersect(as.character(myfiles_SE_all_chr[[32]]$ID), as.character(myfiles_SE_all_chr[[33]]$ID))
int_SE_visual_loss_no_chrXY <- intersect(as.character(myfiles_SE_no_chrXY[[32]]$ID), as.character(myfiles_SE_no_chrXY[[33]]$ID))
int_PE_visual_loss_all_chr <- intersect(as.character(myfiles_PE_all_chr[[32]]$ID), as.character(myfiles_PE_all_chr[[33]]$ID))
int_PE_visual_loss_no_chrXY <- intersect(as.character(myfiles_PE_no_chrXY[[32]]$ID), as.character(myfiles_PE_no_chrXY[[33]]$ID))

# => jaw_claudication                 => 19:21
int_SE_jaw_claudication_all_chr <- intersect(as.character(myfiles_SE_all_chr[[20]]$ID), as.character(myfiles_SE_all_chr[[21]]$ID))
int_SE_jaw_claudication_no_chrXY <- intersect(as.character(myfiles_SE_no_chrXY[[20]]$ID), as.character(myfiles_SE_no_chrXY[[21]]$ID))
int_PE_jaw_claudication_all_chr <- intersect(as.character(myfiles_PE_all_chr[[20]]$ID), as.character(myfiles_PE_all_chr[[21]]$ID))
int_PE_jaw_claudication_no_chrXY <- intersect(as.character(myfiles_PE_no_chrXY[[20]]$ID), as.character(myfiles_PE_no_chrXY[[21]]$ID))

# => ischaemic_features               => 16:18
int_SE_ischaemic_features_all_chr <- intersect(as.character(myfiles_SE_all_chr[[17]]$ID), as.character(myfiles_SE_all_chr[[18]]$ID))
int_SE_ischaemic_features_no_chrXY <- intersect(as.character(myfiles_SE_no_chrXY[[17]]$ID), as.character(myfiles_SE_no_chrXY[[18]]$ID))
int_PE_ischaemic_features_all_chr <- intersect(as.character(myfiles_PE_all_chr[[17]]$ID), as.character(myfiles_PE_all_chr[[18]]$ID))
int_PE_ischaemic_features_no_chrXY <- intersect(as.character(myfiles_PE_no_chrXY[[17]]$ID), as.character(myfiles_PE_no_chrXY[[18]]$ID))

# => GCA_present                      => 1:3
int_SE_GCA_present_all_chr <- intersect(as.character(myfiles_SE_all_chr[[2]]$ID), as.character(myfiles_SE_all_chr[[3]]$ID))
int_SE_GCA_present_no_chrXY <- intersect(as.character(myfiles_SE_no_chrXY[[2]]$ID), as.character(myfiles_SE_no_chrXY[[3]]$ID))
int_PE_GCA_present_all_chr <- intersect(as.character(myfiles_PE_all_chr[[2]]$ID), as.character(myfiles_PE_all_chr[[3]]$ID))
int_PE_GCA_present_no_chrXY <- intersect(as.character(myfiles_PE_no_chrXY[[2]]$ID), as.character(myfiles_PE_no_chrXY[[3]]$ID))

# => Giant_cells                      => 7:9
int_SE_Giant_cells_all_chr <- intersect(as.character(myfiles_SE_all_chr[[8]]$ID), as.character(myfiles_SE_all_chr[[9]]$ID))
int_SE_Giant_cells_no_chrXY <- intersect(as.character(myfiles_SE_no_chrXY[[8]]$ID), as.character(myfiles_SE_no_chrXY[[9]]$ID))
int_PE_Giant_cells_all_chr <- intersect(as.character(myfiles_PE_all_chr[[8]]$ID), as.character(myfiles_PE_all_chr[[9]]$ID))
int_PE_Giant_cells_no_chrXY <- intersect(as.character(myfiles_PE_no_chrXY[[8]]$ID), as.character(myfiles_PE_no_chrXY[[9]]$ID))

# => Media_destruction                => 22:24
int_SE_Media_destruction_all_chr <- intersect(as.character(myfiles_SE_all_chr[[23]]$ID), as.character(myfiles_SE_all_chr[[24]]$ID))
int_SE_Media_destruction_no_chrXY <- intersect(as.character(myfiles_SE_no_chrXY[[23]]$ID), as.character(myfiles_SE_no_chrXY[[24]]$ID))
int_PE_Media_destruction_all_chr <- intersect(as.character(myfiles_PE_all_chr[[23]]$ID), as.character(myfiles_PE_all_chr[[24]]$ID))
int_PE_Media_destruction_no_chrXY <- intersect(as.character(myfiles_PE_no_chrXY[[23]]$ID), as.character(myfiles_PE_no_chrXY[[24]]$ID))

# => Occlusion_grade                  => 28:30
int_SE_Occlusion_grade_all_chr <- intersect(as.character(myfiles_SE_all_chr[[29]]$ID), as.character(myfiles_SE_all_chr[[30]]$ID))
int_SE_Occlusion_grade_no_chrXY <- intersect(as.character(myfiles_SE_no_chrXY[[29]]$ID), as.character(myfiles_SE_no_chrXY[[30]]$ID))
int_PE_Occlusion_grade_all_chr <- intersect(as.character(myfiles_PE_all_chr[[29]]$ID), as.character(myfiles_PE_all_chr[[30]]$ID))
int_PE_Occlusion_grade_no_chrXY <- intersect(as.character(myfiles_PE_no_chrXY[[29]]$ID), as.character(myfiles_PE_no_chrXY[[30]]$ID))

# => Neoangiogenesis                  => 25:27
int_SE_Neoangiogenesis_all_chr <- intersect(as.character(myfiles_SE_all_chr[[26]]$ID), as.character(myfiles_SE_all_chr[[27]]$ID))
int_SE_Neoangiogenesis_no_chrXY <- intersect(as.character(myfiles_SE_no_chrXY[[26]]$ID), as.character(myfiles_SE_no_chrXY[[27]]$ID))
int_PE_Neoangiogenesis_all_chr <- intersect(as.character(myfiles_PE_all_chr[[26]]$ID), as.character(myfiles_PE_all_chr[[27]]$ID))
int_PE_Neoangiogenesis_no_chrXY <- intersect(as.character(myfiles_PE_no_chrXY[[26]]$ID), as.character(myfiles_PE_no_chrXY[[27]]$ID))

# => Intima_pattern                   => 13:15
int_SE_Intima_pattern_all_chr <- intersect(as.character(myfiles_SE_all_chr[[14]]$ID), as.character(myfiles_SE_all_chr[[15]]$ID))
int_SE_Intima_pattern_no_chrXY <- intersect(as.character(myfiles_SE_no_chrXY[[14]]$ID), as.character(myfiles_SE_no_chrXY[[15]]$ID))
int_PE_Intima_pattern_all_chr <- intersect(as.character(myfiles_PE_all_chr[[14]]$ID), as.character(myfiles_PE_all_chr[[15]]$ID))
int_PE_Intima_pattern_no_chrXY <- intersect(as.character(myfiles_PE_no_chrXY[[14]]$ID), as.character(myfiles_PE_no_chrXY[[15]]$ID))

# => Infiltrate_around_vasa_vasorum   => 10:12
int_SE_Infiltrate_around_vasa_vasorum_all_chr <- intersect(as.character(myfiles_SE_all_chr[[11]]$ID), as.character(myfiles_SE_all_chr[[12]]$ID))
int_SE_Infiltrate_around_vasa_vasorum_no_chrXY <- intersect(as.character(myfiles_SE_no_chrXY[[11]]$ID), as.character(myfiles_SE_no_chrXY[[12]]$ID))
int_PE_Infiltrate_around_vasa_vasorum_all_chr <- intersect(as.character(myfiles_PE_all_chr[[11]]$ID), as.character(myfiles_PE_all_chr[[12]]$ID))
int_PE_Infiltrate_around_vasa_vasorum_no_chrXY <- intersect(as.character(myfiles_PE_no_chrXY[[11]]$ID), as.character(myfiles_PE_no_chrXY[[12]]$ID))

##### single-end #####

# Venn diagrams to be made: (fdr)
# => gender (all_chr only)
# => GCA_present (both, all_chr and no_chrXY)
# => Gian_cell (both, all_chr and no_chrXY)
# => Intima_pattern (both)
# => Media_destruction (no_chrXY only => no overlapping)

##### paired-end #####

# Venn diagrams to be made: (fdr)
# => gender (all_chr only)
# => GCA_present (both, all_chr and no_chrXY)
# => Gian_cell (both)
# => Intima_pattern (both)


#venn.diagram(x = list(as.character(myfiles_SE[[32]]$ID), as.character(myfiles_SE[[33]]$ID)), category.names = c("rlog" , "vst"),
#             lwd = 2, lty = 'blank', fill = c("black", "blue"), filename = "visual_loss.png",
#             main = "SE - all_chr - visual_loss")

# save list with transcripts
# "Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES/run_12_Aug20/6_downstream/SE/DESeq2_analysis/transcript_lists"
# "Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/transcript_lists"

# extract string that indicates normalisation method
norm_met_SE_all_chr <- c()
norm_met_SE_no_chrXY <- c()
norm_met_PE_all_chr <- c()
norm_met_PE_no_chrXY <- c()
for (i in 1:length(extracted_names_SE_all_chr)) {
  norm_met_SE_all_chr[i]    <- unlist(str_split(extracted_names_SE_all_chr[[i]][15], "_", n = Inf, simplify = FALSE))[2]
  norm_met_SE_no_chrXY[i]  <- unlist(str_split(extracted_names_SE_no_chrXY[[i]][15], "_", n = Inf, simplify = FALSE))[2]
  norm_met_PE_all_chr[i]    <- unlist(str_split(extracted_names_PE_all_chr[[i]][15], "_", n = Inf, simplify = FALSE))[2]
  norm_met_PE_no_chrXY[i]   <- unlist(str_split(extracted_names_PE_no_chrXY[[i]][15], "_", n = Inf, simplify = FALSE))[2]
}

# create file names
file_names_SE_all_chr <- lapply(extracted_names_SE_all_chr, function(x) paste0(x[10], "_", x[12], "_", x[14]))
file_names_SE_no_chrXY <- lapply(extracted_names_SE_no_chrXY, function(x) paste0(x[10], "_", x[12], "_", x[14]))
file_names_PE_all_chr <- lapply(extracted_names_PE_all_chr, function(x) paste0(x[10], "_", x[12], "_", x[14]))
file_names_PE_no_chrXY <- lapply(extracted_names_PE_no_chrXY, function(x) paste0(x[10], "_", x[12], "_", x[14]))


for (i in 1:length(myfiles_SE_all_chr)){ 
  sink(paste0("/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES/run_12_Aug20/6_downstream/SE/DESeq2_analysis/transcript_lists/transcript_lists_", norm_met_SE_all_chr[i], "_", file_names_SE_all_chr[i], ".txt"))
  cat(as.character(myfiles_SE_all_chr[[i]]$ID), sep = "\n")
  sink()
}

for (i in 1:length(myfiles_SE_no_chrXY)){ 
  sink(paste0("/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES/run_12_Aug20/6_downstream/SE/DESeq2_analysis/transcript_lists/transcript_lists_", norm_met_SE_no_chrXY[i], "_", file_names_SE_no_chrXY[i], ".txt"))
  cat(as.character(myfiles_SE_no_chrXY[[i]]$ID), sep = "\n")
  sink()
}

for (i in 1:length(myfiles_PE_all_chr)){ 
  sink(paste0("/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/transcript_lists/transcript_lists_", norm_met_PE_all_chr[i], "_", file_names_PE_all_chr[i], ".txt"))
  cat(as.character(myfiles_PE_all_chr[[i]]$ID), sep = "\n")
  sink()
}

for (i in 1:length(myfiles_PE_no_chrXY)){ 
  sink(paste0("/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES/run_12_Aug20/6_downstream/PE/DESeq2_analysis/transcript_lists/transcript_lists_", norm_met_PE_no_chrXY[i], "_", file_names_PE_no_chrXY[i], ".txt"))
  cat(as.character(myfiles_PE_no_chrXY[[i]]$ID), sep = "\n")
  sink()
}




# save transformed data and unnormalised as well
#save(dds_all, file="Raw_DESeq_dataset_all.Rda")
#save(vst_all, file="Normalised_DESeq_vst_dataset_all.Rda")
#save(rlog_all, file="Normalised_DESeq_rlog_dataset_all.Rda")

#save(dds_all, file="Raw_DESeq_dataset_noXY.Rda")
#save(vst_all, file="Normalised_DESeq_vst_dataset_noXY.Rda")
#save(rlog_all, file="Normalised_DESeq_rlog_dataset_noXY.Rda")


