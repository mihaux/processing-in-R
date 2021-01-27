# processing script for direct output from Salmon

# source: http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Salmon

# there's also another package that allows to do more processing
# https://bioconductor.org/packages/release/bioc/vignettes/tximeta/inst/doc/tximeta.html

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("tximport", quietly = TRUE)) BiocManager::install("tximport"); library(tximport)

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
output_save <- FALSE

# define directory with data (INPUT)
data_dir <- paste0(main_dir,"/ANALYSES_archived/run_13_Jan21/quants")

# define directory for results (OUTPUT)
dir_out <- paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/processed")
setwd(dir_out)

# get list of .sf all files
files <- list.files(data_dir, full.names = TRUE)

# give proper ID to each sample
names(files) <- paste0("ID_", strtrim(list.files(data_dir,), width = 5))
names(files)[41] <- "ID_8546"

# check if all the input files exist
all(file.exists(files))

# get a list of matrices with original transcript level estimates
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)

# save transcript level counts matrix (175775 x 41)
write.csv(txi.tx$counts, "counts_transcript_level.csv")




# perform gene-level summarization (i.e. from transcript-level to gene-level)
# NOTE: you used an Ensembl transcriptome, so need to use ensembldb packages for gene-level summarization
# http://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html

#txi.sum <- summarizeToGene(txi.tx, tx2gene) 

# save transcript level counts matrix
