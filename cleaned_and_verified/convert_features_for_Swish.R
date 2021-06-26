# convert clinical and histological data spreadsheet from digits to binary yes/no format

# create a new file (to be used for the Swish package) filled with "yes" and "no" only, 
# based on the 'histological_data.csv' file

library(readr)

main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/outliers_excluded"

# load both spreadsheets (without the outliers)
df_hist <- read.csv(file.path(main_dir, "histological_data_no_outliers.csv"), stringsAsFactors = FALSE)
df_clin <- read.csv(file.path(main_dir, "cic_clinical_data_v2_summary_ORDERED_no_outliers.csv"), stringsAsFactors = FALSE)

# exclude sample '16648' whis is NA (might change in the future)
df_hist_new <- df_hist[-which(is.na(df_hist$GCA_present)),]
df_clin_new <- df_clin[-which(df_clin$X=="ID-16648"),]

# EXPLORE THE SPREADSHEETS
# colnames(df_hist_new)
# lapply(df_hist_new, function(x) table(x))

#-----------------------------------------------------------------------------------------#
for (i in 2:ncol(df_hist_new)) {
  cat(colnames(df_hist_new)[i],"\t",table(df_hist_new[,i]), "\n")
}

# Feature name                              "0"     "1"               # at least 5 in each group (for binary)
# GCA_present 	                            7       32                => OK
# Any_granulomatous_infiltrate 	            37      2                 => to be dropped
# Granulomatous_infiltrate_in_adventitia 	  37      2                 => to be dropped
# Granulomatous_infiltrate_in_media 	      39      0                 => to be dropped
# Granulomatous_infiltrate_in_intima 	      38      1                 => to be dropped 
# Any_lymphocytic_infiltrate 	              7       32                => OK 
# Lymphocytic_infiltrate_in_adventitia 	    7       32                => OK 
# Lymphocytic_infiltrate_in_media 	        12      27                => OK 
# Lymphocytic_infiltrate_in_intima 	        8       31                => OK 
# Barcelona_score 	                        7 1 6 25    (multi-class) => to be dropped
# Adventitia_pattern 	                      7 1 7 24    (multi-class) => OK
# Media_pattern 	                          14 7 12 6   (multi-class) => OK
# Intima_pattern 	                          8 8 15 8    (multi-class) => OK
# Giant_cells 	                            29      10                => OK 
# Aggregates 	                              35      4                 => to be dropped 
# Infiltrate_around_vasa_vasorum 	          22      17                => OK 
# PALI 	                                    35      4                 => to be dropped
# Media_destruction 	                      19      20                => OK 
# Neoangiogenesis 	                        28      11                => OK 
# Hyperplasia 	                            7       32                => OK 
# Fibrosis 	                                8       31                => OK 
# Oedema 	                                  35      4                 => to be dropped 
# Occlusion_grade 	                        7 4 6 21 1  (multi-class) => OK

# exclude the sample with NA in GCA_present
df_hist_fin <- df_hist_new[-which(is.na(df_hist$GCA_present)),]

# subset both dataframes to retain only features of interest
df_hist_fin <- df_hist_new[,c("Database_number", "GCA_present", "Any_lymphocytic_infiltrate", 
                              "Lymphocytic_infiltrate_in_adventitia","Lymphocytic_infiltrate_in_media", 
                              "Lymphocytic_infiltrate_in_intima", "Adventitia_pattern", "Media_pattern", 
                              "Intima_pattern", "Giant_cells", "Infiltrate_around_vasa_vasorum", 
                              "Media_destruction", "Neoangiogenesis", "Hyperplasia", "Fibrosis", "Occlusion_grade")]

#-----------------------------------------------------------------------------------------#
for (i in 8:ncol(df_clin_new)) {
  cat(colnames(df_clin_new)[i],"\t",table(df_clin_new[,i]), "\n")
}

# NOTE: only selected features from the output kept below
# Batch..Ian. 	                            25      14                => to be dropped (for now)
# visual_loss 	                            24      15                => OK
# jaw_claudication 	                        15      24                => OK 
# ischaemic_features 	                      10      29                => OK
# gender 	                                  15      24                => OK
# Multi-class features
# year.TAB.sample.was.collected 	 2 7 18 11 1                        => to be dropped (for now)
# number.of.days.on.steroids.at.TAB 	 1 2 1 7 5 3 6 8 2 1 1 2        => OK
# number.of.days.between.TAB.and.BL.blood.sample 	 1 5 7 1 1 5 2 4 3 6 2 1 1 => to be dropped (for now)
# age.at.BL 	 1 1 1 4 2 2 2 2 3 1 5 1 3 1 3 1 1 2 1 1 1              => OK

df_clin_fin <- df_clin_new[,c("X", "visual_loss", "jaw_claudication", "ischaemic_features", 
                              "gender", "number.of.days.on.steroids.at.TAB", "age.at.BL")]
#-----------------------------------------------------------------------------------------#

# check if IDs match between df_hist_fin and df_clin_fin
df_hist_fin$Database_number == df_clin_fin$X

#---------------------------------------------------------------------#
# create 3 functions that will replace with 0 and 1 (could merge them all in one function)

# convert_0_1() function works for all binary features
# => 0 => replace with "no" 
# => 1 => replace with "yes"
# ARGS: (1) df - data.frame, (2) list_of_feats - feature name (also column name)
# NOTE: list_of_feats needs to include the column with row names at first position

convert_0_1 <- function(df, list_of_feats){
  temp_0 <- df[,list_of_feats]
  temp_1 <- lapply(temp_0, function(x) replace(x, x==0, "no"))
  temp_2 <- lapply(temp_1, function(x) replace(x, x==1, "yes"))
  return(as.data.frame(temp_2))
}

# convert_0_3() function works for: Adventitia_pattern, Media_pattern and Intima_pattern
# => 0 and 1 => replace with "no" 
# => 2 and 3 => replace with "yes"
# ARGS: (1) df - data.frame, (2) feat - feature name (also column name)

convert_0_3 <- function(df, list_of_feats){
  temp_0 <- df[,list_of_feats]
  temp_1 <- lapply(temp_0, function(x) replace(x, x==0, "no"))
  temp_2 <- lapply(temp_1, function(x) replace(x, x==1, "no"))
  temp_3 <- lapply(temp_2, function(x) replace(x, x==2, "yes"))
  temp_4 <- lapply(temp_3, function(x) replace(x, x==3, "yes"))
  return(as.data.frame(temp_4))
}

# convert_0_4() function works for: Occlusion_grade and Barcelona_score (DO NOT USE FOR NOW)
# => less than 3            => replace with 0
# => equal or more than 3   => replace with 1
convert_0_4 <- function(df, list_of_feats){
  temp_0 <- df[,list_of_feats]
  temp_1 <- lapply(temp_0, function(x) replace(x, x==0, "no"))
  temp_2 <- lapply(temp_1, function(x) replace(x, x==1, "no"))
  temp_3 <- lapply(temp_2, function(x) replace(x, x==2, "no"))
  temp_4 <- lapply(temp_3, function(x) replace(x, x==3, "yes"))
  temp_5 <- lapply(temp_4, function(x) replace(x, x==4, "yes"))
  return(as.data.frame(temp_5))
}

convert_gender <- function(df, list_of_feats){
  temp_0 <- df[,list_of_feats]
  temp_1 <- lapply(temp_0, function(x) replace(x, x==1, "no"))    # male
  temp_2 <- lapply(temp_1, function(x) replace(x, x==2, "yes"))   # female
  return(as.data.frame(temp_2))
}

# NOTE: the ">=" operator applies to characters as well, so it needs to be used first

# age.at.BL   < 75 y/o    |   ≥ 75 y/o
convert_age <- function(df, feat){    # feat=age.at.BL
  temp_0 <- df[,feat]
  temp_1 <- replace(temp_0, temp_0>=75, "yes")   # ≥ 75 y/o
  temp_2 <- replace(temp_1, temp_1<75, "no")    # < 75 y/o
  return(as.data.frame(temp_2))
}

# number.of.days.on.steroids.at.TAB  < 6 days    |   ≥ 6 days
convert_steroids <- function(df, feat){ # feat=number.of.days.on.steroids.at.TAB
  temp_0 <- df[,feat]
  temp_1 <- replace(temp_0, temp_0>=6, "yes")   # ≥ 6 days
  temp_2 <- replace(temp_1, temp_1<6, "no")     # < 6 days
  return(as.data.frame(temp_2))
}

### histological ###
# binary: 
# GCA_present | Giant_cells | Infiltrate_around_vasa_vasorum | Media_destruction | Neoangiogenesis

# multi-class:
# Adventitia_pattern | Media_pattern | Intima_pattern => feature | 0: normal  | 1: focal  | 2: multifocal | 3: diffuse
# Occlusion_grade                                     => feature | 0: none    | 1: <50%   | 2: 50-75%     | 3: 75-100%  | 4: complete

# for binary only
df_hist_final_I <- convert_0_1(df_hist_fin, c("Database_number", "GCA_present", "Any_lymphocytic_infiltrate", 
                                              "Lymphocytic_infiltrate_in_adventitia","Lymphocytic_infiltrate_in_media", 
                                              "Lymphocytic_infiltrate_in_intima", "Giant_cells", "Infiltrate_around_vasa_vasorum", 
                                              "Media_destruction", "Neoangiogenesis", "Hyperplasia", "Fibrosis"))

# for multi-class only
df_hist_final_II <- convert_0_3(df_hist_fin, c("Database_number", "Adventitia_pattern", "Media_pattern", "Intima_pattern"))
df_hist_final_III <- convert_0_4(df_hist_fin, c("Database_number", "Occlusion_grade"))

# merge all dataframes for histological features
df_hist_final_all <- data.frame(df_hist_final_III, df_hist_final_I[,-1], df_hist_final_II[,-1])
names(df_hist_final_all)[1] <- "names"

### clinical ###
# binary: visual_loss | jaw_claudication | ischaemic_features | gender (special case)
# for binary only
df_clin_final_I <- convert_0_1(df_clin_fin, c("X", "visual_loss", "jaw_claudication", "ischaemic_features"))

# special case: gender    1 (male) | 2 (female)
df_clin_final_II <- convert_gender(df_clin_fin, c("X", "gender"))

# continuous:
# number of days on steroids  < 6 days    |   ≥ 6 days
# age                         < 75 y/o    |   ≥ 75 y/o

df_clin_final_III <- convert_age(df_clin_fin, c("age.at.BL"))
df_clin_final_IV <- convert_steroids(df_clin_fin, c("number.of.days.on.steroids.at.TAB"))

# merge all dataframes for clinical features
df_clin_final_all <- data.frame(df_clin_final_I, df_clin_final_II[-1], df_clin_final_III, df_clin_final_IV)
names(df_clin_final_all) <- c("names", "visual_loss", "jaw_claudication", "ischaemic_features",
                              "gender", "age", "steroids_duration")
  
# FINAL OUTPUT FILES:
# df_clin_final_all
# df_hist_final_all

# merge both dataframes
df_final_both <- data.frame(df_clin_final_all, df_hist_final_all[-1])

# write all output files
out_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/sample_info_for_swish_no_outliers"
setwd(out_dir)

write_delim(df_clin_final_all, "samples_clinical.txt", delim="\t")
write_delim(df_hist_final_all, "samples_histological.txt", delim="\t")
write_delim(df_final_both, "samples_all.txt", delim="\t")

# get number of YES and NO for each feature
summary_yes_no <- t(apply(df_final_both[,-1], 2, table))

summary_yes_no_final <- cbind(rownames(summary_yes_no), summary_yes_no)
colnames(summary_yes_no_final) <- c("Feature_name", "no", "yes")

write_delim(as.data.frame(summary_yes_no_final), "summary_yes_no.txt", delim="\t")



