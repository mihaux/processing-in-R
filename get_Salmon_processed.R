# processing script for direct output from Salmon

# source: http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Salmon

# on twitter: https://twitter.com/mikelove/status/1163783031741702144 
# on github: https://github.com/mikelove/swish-demos/blob/master/differential-transcript-airway.knit.md

# https://rdrr.io/bioc/fishpond/man/fishpond-package.html

# there's also another package that allows to do more processing (hard to get it installed)
# https://bioconductor.org/packages/release/bioc/vignettes/tximeta/inst/doc/tximeta.html

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("tximport", quietly = TRUE)) BiocManager::install("tximport"); library(tximport)      # for transcript-level
if (!requireNamespace("ensembldb", quietly = TRUE)) BiocManager::install("ensembldb"); library(ensembldb)   # for gene-level

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
data_dir <- paste0(main_dir,"/ANALYSES_archived/run_13_Jan21/quant_new_idx")

# define directory for results (OUTPUT)
dir_out <- paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/processed_inf_rep")
setwd(dir_out)

# get list of .sf all files
files <- list.files(data_dir, full.names = TRUE)

# give proper ID to each sample
names(files) <- paste0("ID_", strtrim(list.files(data_dir,), width = 5))
names(files)[41] <- "ID_8546"

# check if all the input files exist
all(file.exists(files))

########## transcript-level ##########


# get a list of matrices with original transcript-level estimates
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)

txi.inf.rep <- tximport(files, type = "salmon", txOut = TRUE)
names(txi.inf.rep)
head(txi.inf.rep$counts)

txi <- tximport(files, type = "none", txOut = TRUE, txIdCol = "Name", abundanceCol = "TPM", 
                countsCol = "NumReads", lengthCol = "Length", importer = function(x) read_tsv(x, skip = 8))

# generate counts from abundances, scaled to library size, "scaledTPM"
txi.tx_scaledTPM <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = c("scaledTPM"), dropInfReps = FALSE)

# generate counts from abundances, scaled using the average transcript length, averaged over samples and to library size, "lengthScaledTPM"
txi.tx_lengthScaledTPM <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"))


# save transcript-level counts matrix (175775 x 41)
if(output_save==TRUE){
  write.csv(txi.tx$counts, "counts_transcript-level.csv")
  write.csv(txi.tx_scaledTPM$counts, "counts_transcript-level_scaledTPM.csv")
  write.csv(txi.tx_lengthScaledTPM$counts, "counts_transcript-level_lengthScaledTPM.csv")
}

# exclude the outliers ("ID_8546") from transcript-level counts matrix
outliers_trans <- which(colnames(txi.tx$counts)=="ID_8546")

if(output_save==TRUE){
  write.csv(txi.tx$counts[,-outliers_trans], "counts_transcript-level_no_outliers.csv")
  write.csv(txi.tx_scaledTPM$counts[,-outliers_trans], "counts_transcript-level_scaledTPM_no_outliers.csv")
  write.csv(txi.tx_lengthScaledTPM$counts[,-outliers_trans], "counts_transcript-level_lengthScaledTPM_no_outliers.csv")
}

########## gene-level ##########
# NOTE: you used an Ensembl transcriptome, so need to use ensembldb packages for gene-level summarization
# http://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html

# load annotation package (make sure the version is correct)
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) BiocManager::install("EnsDb.Hsapiens.v86")   # GRCh37 (latest)
#if (!requireNamespace("EnsDb.Hsapiens.v79", quietly = TRUE)) BiocManager::install("EnsDb.Hsapiens.v79")  # GRCh38 
#if (!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE)) BiocManager::install("EnsDb.Hsapiens.v75")  # GRCh37

# NOTE: Transcripts need to be associated with gene IDs for gene-level summarization.

# create a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. 

library(EnsDb.Hsapiens.v86)
txdb <- EnsDb.Hsapiens.v86

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# NOTE: the tird column ("TXID") is redundant
tx2gene$TXID <- NULL

# perform gene-level summarization (i.e. from transcript-level to gene-level)
txi.sum <- summarizeToGene(txi.tx, tx2gene, ignoreTxVersion=TRUE) # by default the argument countsFromAbundance = c("no")

# in txi.sum (a list of 4 matrices)
# "abundance" - provided by the quantification tools as TPM (transcripts-per-million)
# "counts"    - estimated counts (possibly fractional)
# "length"    - contains the effective gene lengths
# NOTE: The "length" matrix can be used to generate an offset matrix for downstream gene-level differential analysis of count matrices.

# generate counts from abundances, scaled to library size, "scaledTPM"
txi.sum_scaledTPM <- summarizeToGene(txi.tx, tx2gene, ignoreTxVersion=TRUE, countsFromAbundance = c("scaledTPM")) 

# generate counts from abundances, scaled using the average transcript length, averaged over samples and to library size, "lengthScaledTPM"
txi.sum_lengthScaledTPM <- summarizeToGene(txi.tx, tx2gene, ignoreTxVersion=TRUE, countsFromAbundance = c("lengthScaledTPM")) 

# NOTE: Using either of these approaches, the counts are not correlated with length, and so the length matrix should not be provided as an offset for downstream analysis packages. 

# argument: countsFromAbundance
# Whether to generate estimated counts using abundance estimates:
# => scaled up to library size (scaledTPM),
# => scaled using the average transcript length over samples and then the library size (lengthScaledTPM), or
# => scaled using the median transcript length among isoforms of a gene, and then the library size (dtuScaledTPM).

# dtuScaledTPM doesn't work
# NOTE: dtuScaledTPM is designed for DTU analysis in combination with txOut=TRUE, and it requires specifing a tx2gene data.frame. 
# dtuScaledTPM works such that within a gene, values from all samples and all transcripts get scaled by the same fixed median transcript length. 
# If using scaledTPM, lengthScaledTPM, or geneLengthScaledTPM, the counts are no longer correlated across samples with transcript length, and so the length offset matrix should not be used.

# save gene-level counts matrices (37 361 x 41)
if(output_save==TRUE){
  write.csv(txi.sum$counts, "counts_gene-level.csv")
  write.csv(txi.sum_scaledTPM$counts, "counts_gene-level_scaledTPM.csv")
  write.csv(txi.sum_lengthScaledTPM$counts, "counts_gene-level_lengthScaledTPM.csv")
}

# exclude the outliers ("ID_8546") from gene-level counts matrix
outliers_gene <- which(colnames(txi.sum$counts)=="ID_8546")

# save gene-level counts matrices (37 361 x 40)
if(output_save==TRUE){
  write.csv(txi.sum$counts[,-outliers_gene], "counts_gene-level_no_outliers.csv")
  write.csv(txi.sum_scaledTPM$counts[,-outliers_gene], "counts_gene-level_scaledTPM_no_outliers.csv")
  write.csv(txi.sum_lengthScaledTPM$counts[,-outliers_gene], "counts_gene-level_lengthScaledTPM_no_outliers.csv")
}




library(tximportData)
dir <- system.file("extdata", package = "tximportData")
list.files(dir)

samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)

files <- file.path(dir, "salmon_gibbs", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
txi.inf.rep <- tximport(files, type = "salmon", txOut = TRUE)
names(txi.inf.rep)







