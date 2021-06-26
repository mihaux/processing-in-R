# this script processes slide scores spreadsheet by selecting the highest score across all sections for each feature and each sample

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(readr)

# INPUT: /Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/slide_scores_v5.csv
# OUTPUT: /Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/slide_scores_v6.csv

#args <- commandArgs(trailingOnly = TRUE)

args <- c("/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/histological_data/Image_atlas_data_MODIFIED_extracted.csv",
          "/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/histological_data")

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

# load slide scores table
df <- read.csv(args[1], header = TRUE)
  

### NOT IN USE ANYMORE ###
# empty columns are replaced with NA (i.e. these columns where GCA.present. = 0 or "NS")
# replace NA with zeros 
#which(df_subset$GCA.present. == 0)      # 26 rows
#which(df_subset$GCA.present. == "NS")   # 13 rows

#for(i in 1:length(df_subset$GCA.present.)){
#  if(df_subset$GCA.present.[i] == 0 || df_subset$GCA.present.[i] == "NS"){
#    df_subset[i,3:ncol(df_subset)] <- 0
#  }
#}

# there's one exeption that needs to be handled manually
#df_subset[which(df_subset$GCA.present. == "0 see PALI"),]
# this should have 0 in GCA present and 1 in the PALI column – please retain data from section 3
# => keep column 19 as 1 and the rest as 0s
# set 0 from col 3 to 18 and 20 to 25
#which(df_subset$GCA.present. == "0 see PALI") # => [1] 6
#df_subset[6,3:18] <- 0
#df_subset[6,20:25] <- 0

# there's also another missing value: ID - 8546; section - 3 => df_subset[3,]; change to 0
#df_subset[3,25] <- 0

##########################

# get indices of rows for each samples
all_subsets <- lapply(unique(df$Database_number), function(x) which(df$Database_number == x))
names(all_subsets) <- unique(df$Database_number)

# divide df_subset into 41 small dataframes including only one sample
all_subsets_final <- list()
for(i in 1:length(all_subsets)){
  all_subsets_final[[i]] <- df[all_subsets[[i]],]
}

# for each of these 41 dataframes, create a function that will read each df by column and add a new row with the highest score in current column

# create a testing dataframe
#df_test <- all_subsets_final[1:3]

outputs_vec <- c()
outputs_lis <- list()
# read all_subsets_final list one by one (one small dataframe at a time)
for(i in 1:length(all_subsets_final)){
  for(j in 1:ncol(all_subsets_final[[1]])){
    outputs_vec[j] <- max(as.vector(all_subsets_final[[i]][,j]), na.rm = TRUE)
    print(max(as.vector(all_subsets_final[[i]][,j])))
    
  }
  outputs_lis[[i]] <- outputs_vec
}

# convert list to dataframe and add col & row names
df_final <- data.frame(Reduce(rbind, outputs_lis))
colnames(df_final) <- colnames(df)
rownames(df_final) <- 1:41

# add "highest score" in each row in 'Section.number' column
#df_final$Section_number <- rep("highest score", 41)
df_final$Section_number <- NULL

# write final dataframe
write_delim(df_final, "histological_data_bis.csv", delim = ",")


