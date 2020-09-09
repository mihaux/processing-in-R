# this script scans metadata tables and return the number of patients per group for each feature
# NOTE: run separately for each .csv file

# INPUT:  .csv tables with clinical and histological features
# OUTPUT: .txt file with a report

# get working directory to recognise the machine
w_dir <- getwd()

# create a shortcut for the OneDrive directory where all files are stored
if(startsWith(w_dir, "/Users/michal")){           
  main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"    # on my mac
  cat("Running on Michal's MacBookPro")
} else if (startsWith(w_dir, "/Users/ummz")) {    
  main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds"      # on uni mac  
  cat("Running on med-mac-211856")
} else {
  print("Unrecognised machine.")
}

# args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input: clinical) path to .csv with clinical data,
       \n(2 - input: histological) path to .csv with histological data and,
       \n(3 - output) path where output file should be stored", call.=FALSE)
}

args <- c(paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/data/metadata/slide_scores/slide_scores_v6.csv"),
          paste0(main_dir, "/ANALYSES/run_12_Aug20/6_downstream/"))

cat("Directory with metadata (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load data 
df_clinical <- read.csv(args[1], row.names = 1, header = TRUE)        # dim = 41 19
df_histological <- read.csv(args[2], row.names = 1, header = TRUE)    # dim = 40 24

if(FALSE){ # there are 4 sheets in cic_clinical_data_v2.csv | the other 3 are not important for now
  # load each sheet from cic_clinical_data_v2.xlsx
  setwd("/Users/ummz/Documents/OneDrive - University of Leeds/data/metadata/cic_clinical_data_v2_split/")
  df_clinical_sheet_1 <- read.csv("cic_clinical_data_v2_summary.csv")                               # 41 x 19
  df_clinical_sheet_2 <- read.csv("cic_clinical_data_v2_TABUL_clinical_data.csv", header = FALSE)   # 36 x 142; there's double header
  df_clinical_sheet_3 <- read.csv("cic_clinical_data_v2_TABUL_steroid_data.csv")                    # 67 x 6
  df_clinical_sheet_4 <- read.csv("cic_clinical_data_v2_UK_GCA_clinical_data.csv")                  # 7 x 149 => unimportant 
  
  # OTHER ISSUES:
  
  # there's double header so it need to be fixed
  #sheet_2_mod <- sheet_2[-c(1,2),]
  
  #mat_temp <- as.matrix(sheet_2)
  #colnames(sheet_2_mod) <- as.character(mat_temp[2,])
  
  # match 'study ID' with 'database code' from Summary
  # keep only samples that are present in sheet_1 (i.e. from my cohort)
  # it includes 35 common samples with sheet_1
  # 6 are missing: 12331, 12330, 14455, 11026, 11028, 8546
  
  # add "database.code" that matches with `Study ID`

  # match 'study ID' with 'database code' from Summary
  # keep only samples that are present in sheet_1 (i.e. from my cohort)
  #to_be_kept_s3 <- which(as.character(sheet_1$Study.ID) %in% as.character(sheet_3$Study.ID)) # 20
  
  #sheet_3_final <- sheet_3[to_be_kept_s3,]
  
  # there are duplicated rows; all study ID's in sheet_3 are included in sheet_1
  #unique(sheet_3$Study.ID) == intersect(as.character(sheet_1$Study.ID), as.character(sheet_3$Study.ID))
}

# add "ID_" to all rownames
rownames(df_clinical) <- paste0("ID_", rownames(df_clinical))
rownames(df_histological) <- paste0("ID_", rownames(df_histological))

cat("Total number of patients in clinical data:", nrow(df_clinical), "\n")
cat("Total number of patients in histological data:", nrow(df_histological), "\n")

cat("Missing sample in the histological data:",  setdiff(rownames(df_clinical), rownames(df_histological)), "\n")

# print all features in each spreadsheet
cat("All clinical features:\n")
colnames(df_clinical)

cat("All histological features:\n")
colnames(df_histological)

### Features in cic_clinical_data_v2.xlsx (summary sheet) ###

## (1) check how many class (levels) for each feature
out_clinical <- list()
for (i in 1:ncol(df_clinical)) {
  out_clinical[i] <- unique(df_clinical[i])
}

#unlist(lapply(out_clinical, length))

# columns 1-7 can be ignored as they don't contain relevant information

# create new objects where 1-7 columns will not be included.
df_clinical_bis <- df_clinical[-c(1:7)]
out_clinical_bis <- out_clinical[-c(1:7)]

## (2) check how many cases in each class?
for (i in 1:ncol(df_clinical_bis)) {
  cat("###", colnames(df_clinical_bis[i]), "###\n")
  for (j in 1:length(out_clinical_bis[[i]])) {
    cat(out_clinical_bis[[i]][j], "=> ", length(which(as.character(out_clinical_bis[[i]][j]) == df_clinical_bis[i])), "\n")
  }
  cat("\n")
}

if(FALSE){
  # Visual loss at BL = reduced or lost vision (temporary or permanent) pre steroids or at BL assessment or AION or PION 
  # Jaw claudication at BL = present pre steroids or at baseline assessment visit
  # Ischaemic features at BL = jaw claudication, tongue claudication or visual loss (temporary or permanent) pre steroids or at baseline
  
  # for columns 3-9: 0 - NO | 1 - YES
  # for columns 5-6: 99 - not assessed
  # for column 9: 2 - unknown
  # for column 10: 1 - male | 2 - female
}

### Features in slide_scores_v6.csv (summary sheet) ###

## (1) check how many class (levels) for each feature
out_histological <- list()
for (i in 1:ncol(df_histological)) {
  out_histological[i] <- unique(df_histological[i])
}

#unlist(lapply(out_histological, length))

## (2) check how many cases in each class?
for (i in 1:ncol(df_histological)) {
  cat("###", colnames(df_histological[i]), "###\n")
  for (j in 1:length(out_histological[[i]])) {
    cat(out_histological[[i]][j], "=> ", length(which(as.character(out_histological[[i]][j]) == df_histological[i])), "\n")
  }
  cat("\n")
}
