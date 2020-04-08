# this script processes slide scores spreadsheet by selcting ony those that are in my cohort and then 

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# INPUT: /Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/slide_scores_v5.csv
# OUTPUT: /Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/slide_scores_v6.csv

#args <- commandArgs(trailingOnly = TRUE)

args <- c("/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/slide_scores_v5.csv",
          "/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/")

cat("Example of usage: \n Rscript metadata_proc.R /Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/slide_scores_v5.csv /Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/")

if (length(args)!=2) {
  stop("2 arguments must be supplied: 
       \n(1 - input) path to .csv file with metadata (slide scores), 
       \nand (2 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[1], sep="\n")

cat("Directory for results (OUT): ")
cat(args[2], sep="\n")

setwd(args[2])

# TODO:
# add HS (highest documented score) for each row of the resulting table in 'Section.number' column 
# modify column names 12, 13, 14, 15 and 25 (see below)
# if the score for GCA.present is zero then put a 0 in every column
# ATTENTION: there is one exception which is 12331 – this should have 0 in GCA present and 1 in the PALI column – please retain data from section 3
# 16648 is the one with no scores. 
# There are two entries for 8546 and 14914 and possibly others.
# sort on database number 

# COLUMN NAMES
# [1] "Database.number"                                                                                                             
# [2] "Section.number"                                                                                                  
# [3] "GCA.present."                                                                                                
# [4] "Any.granulomatous.infiltrate."                                                                               
# [5] "Granulomatous.infiltrate.in.adventitia."                                                                     
# [6] "Granulomatous.infiltrate.in.media."                                                                          
# [7] "Granulomatous.infiltrate.in.intima."                                                                         
# [8] "Any.lymphocytic.infiltrate."                                                                                 
# [9] "X"                                                                                                           
# [10] "Lymphocytic.infiltrate.in.media."                                                                            
# [11] "Lymphocytic.infiltrate.in.intima."                                                                           
# [12] "Barcelona.score."                       => DROP: ..0..normal..1..adventitial..2..adventitial.invasive..3..concentric.bilayer..4..panarteritis."
# [13] "Adventitia.pattern."                    => DROP: ..0..normal..1..focal..2..multifocal..3..diffuse."                                         
# [14] "Media.pattern."                         => DROP: ..0..normal..1..focal..2..multifocal..3..diffuse."                                              
# [15] "Intima.pattern."                        => DROP: ..0..normal..1..focal..2..multifocal..3..diffuse."                                             
# [16] "Giant.cells."                                                                                                
# [17] "Aggregates."                                                                                                 
# [18] "Infiltrate.around.vasa.vasorum."                                                                             
# [19] "PALI."                                                                                                       
# [20] "Media.destruction."                                                                                          
# [21] "Neoangiogenesis."                                                                                            
# [22] "Hyperplasia."                                                                                                
# [23] "Fibrosis."                                                                                                   
# [24] "Oedema."                                                                                                     
# [25] "Occlusion.grade."                       => DROP: ..0..none..1...50...2..50.75...3..75.100...4..complete."  

# load slide scores table
df <- read.csv(args[1], header = TRUE)

# NOTE: colnames(df)[9] = missing !!!
# change colnames where needed
colnames(df)[12] <- "Barcelona.score."
colnames(df)[13] <- "Adventitia.pattern."
colnames(df)[14] <- "Media.pattern."
colnames(df)[15] <- "Intima.pattern."
colnames(df)[25] <- "Occlusion.grade."

# read IDsas vector
IDs <- scan('sample_IDs.txt', character(), quote = "")

# there are some "Database.number" that are incorrect
# TODO: check if they are included in you cohort: 
# 14714, 14767, 14883, 14933, 14581, 14582, 14583, 14584, 14585, 14586, 14587, 14588, 14708, 14709, 14H15502

# 176 => [1] 14714  7-9
# 179 => [1] 14714  4-6
# 182 => [1] 14714 1-3
# 184 => [1] 14767 7-9
# 186 - 188 => [1] 14767 4-6 
#                  14767 1-3 
#                  14883 7-9
# 190 - 192 => [1] 14883 4-6 
#                  14883 1-3 
#                  14933 7-9
# 195 - 212 => [1] 14933 4-6 
#                  14933 1-3 
#                  14581 1-3 
#                  14582 1-3 
#                  14583 1-3 
#                  14582 1-3 
#                  14583 1-3 
#                  14584 1-3 
#                  14585 1-3
#                  14586 1-3 
#                  14583 1-3 
#                  14584 1-3 
#                  14585 1-3 
#                  14586 1-3 
#                  14587 1-3 
#                  14588 1-3 
#                  14708 1-3 
#                  14709 1-3
# 378 => [1] 14H15502

to_be_checked <- c(14714, 14767, 14883, 14933, 14581, 14582, 14583, 14584, 14585, 14586, 14587, 14588, 14708, 14709)
additional <- "14H15502" # => Wrong tissue

intersect(IDs, to_be_checked) # => numeric(0) 
# CONCLUSION: none of these incorrect IDs are included in my cohort

setdiff(IDs, df$Database.number) # => [1] "14058" => this sample is not included in slide scores spreadsheet

IDs_updated <- IDs[-which(IDs == 14058)]
  
# get indices of rows to be kept (i.e. those from my cohort)
to_be_extracted <- which(df$Database.number %in% IDs_updated)

# create a subset including only samples from your cohort
df_subset <- df[to_be_extracted,]

# save df_subset
# write.csv(df_subset, "subset_temp_1.csv")

# empty columns are replaced with NA (i.e. these columns where GCA.present. = 0 or "NS")
# replace NA with zeros 
which(df_subset$GCA.present. == 0)      # 26 rows
which(df_subset$GCA.present. == "NS")   # 13 rows

for(i in 1:length(df_subset$GCA.present.)){
  if(df_subset$GCA.present.[i] == 0 || df_subset$GCA.present.[i] == "NS"){
    df_subset[i,3:ncol(df_subset)] <- 0
  }
}

# there's one exeption that needs to be handled manually
df_subset[which(df_subset$GCA.present. == "0 see PALI"),]
# this should have 0 in GCA present and 1 in the PALI column – please retain data from section 3
# => keep column 19 as 1 and the rest as 0s
# set 0 from col 3 to 18 and 20 to 25
which(df_subset$GCA.present. == "0 see PALI") # => [1] 6
df_subset[6,3:18] <- 0
df_subset[6,20:25] <- 0

# save df_subset
# write.csv(df_subset, "subset_temp_2.csv")

# there's also another missing value: ID - 8546; section - 3 => df_subset[3,]; change to 0
df_subset[3,25] <- 0

# get indices of rows for each samples
all_subsets <- lapply(IDs_updated, function(x) which(df_subset$Database.number == x))
names(all_subsets) <- IDs_updated

# divide df_subset into 41 small dataframes including only one sample
all_subsets_final <- list()
for(i in 1:length(all_subsets)){
  all_subsets_final[[i]] <- df_subset[all_subsets[[i]],]
}

# for each of these 41 dataframes, create a function that will read each df by column and add a new row with the highest score in current column

# create a testing dataframe
#df_test <- all_subsets_final[1:3]

outputs_vec <- c()
outputs_lis <- list()
# read all_subsets_final list one by one (one small dataframe at a time)
for(i in 1:length(all_subsets_final)){
  for(j in 1:ncol(all_subsets_final[[1]])){
    outputs_vec[j] <- max(as.vector(all_subsets_final[[i]][,j]))
  }
  outputs_lis[[i]] <- outputs_vec
}

# convert list to dataframe and add col & row names
df_final <- data.frame(Reduce(rbind, outputs_lis))
colnames(df_final) <- colnames(df_subset)
rownames(df_final) <- 1:40

# add "highest score" in each row in 'Section.number' column
df_final$Section.number <- rep("highest score", 40)

# write final dataframe
write.csv(df_final, "subset_final.csv", row.names = FALSE)

