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
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2"); suppressMessages(library(DESeq2))

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
data_dir <- paste0(main_dir,"/ANALYSES_archived/run_13_Jan21/quants")

# define directory for results (OUTPUT)
dir_out <- paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/")
setwd(dir_out)

# get list of .sf all files
files <- list.files(data_dir, full.names = TRUE)

# give proper ID to each sample
names(files) <- paste0("ID-", strtrim(list.files(data_dir,), width = 5))
names(files)[41] <- "ID-8546"

# check if all the input files exist
all(file.exists(files))

########## transcript-level ##########
# get a list of matrices with original transcript-level estimates
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)

# generate counts from abundances, scaled to library size, "scaledTPM"
txi.tx_scaledTPM <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = c("scaledTPM"))

# generate counts from abundances, scaled using the average transcript length, averaged over samples and to library size, "lengthScaledTPM"
txi.tx_lengthScaledTPM <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"))

# save transcript-level counts matrix (175775 x 41)
if(output_save==TRUE){
  write.csv(txi.tx$counts, "processed/counts_transcript-level.csv", )
  write.csv(txi.tx_scaledTPM$counts, "processed/counts_transcript-level_scaledTPM.csv")
  write.csv(txi.tx_lengthScaledTPM$counts, "processed/counts_transcript-level_lengthScaledTPM.csv")
}

# exclude the outliers ("ID_8546") from transcript-level counts matrix
outliers_trans <- which(colnames(txi.tx$counts)=="ID-8546")

if(output_save==TRUE){
  write.csv(txi.tx$counts[,-outliers_trans], "processed/counts_transcript-level_no_outliers.csv")
  write.csv(txi.tx_scaledTPM$counts[,-outliers_trans], "processed/counts_transcript-level_scaledTPM_no_outliers.csv")
  write.csv(txi.tx_lengthScaledTPM$counts[,-outliers_trans], "processed/counts_transcript-level_lengthScaledTPM_no_outliers.csv")
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
  write.csv(txi.sum$counts, "processed/counts_gene-level.csv")
  write.csv(txi.sum_scaledTPM$counts, "processed/counts_gene-level_scaledTPM.csv")
  write.csv(txi.sum_lengthScaledTPM$counts, "processed/counts_gene-level_lengthScaledTPM.csv")
}

# exclude the outliers ("ID_8546") from gene-level counts matrix
outliers_gene <- which(colnames(txi.sum$counts)=="ID-8546")

# save gene-level counts matrices (37 361 x 40)
if(output_save==TRUE){
  write.csv(txi.sum$counts[,-outliers_gene], "processed/counts_gene-level_no_outliers.csv")
  write.csv(txi.sum_scaledTPM$counts[,-outliers_gene], "processed/counts_gene-level_scaledTPM_no_outliers.csv")
  write.csv(txi.sum_lengthScaledTPM$counts[,-outliers_gene], "processed/counts_gene-level_lengthScaledTPM_no_outliers.csv")
}

# N O R M A L I S A T I O N

# define directory for results selected for use in statisitcal testing
dir_out_norm <- paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/statistical_testing/input")
setwd(dir_out)

# NOTE: use "scaledTPM" counts as raw and perform some sort of normlisation (using edgeR package)

### =>     RAW     <= ###
cts_raw_transcript <- txi.tx_scaledTPM$counts     # 175 775 x 41  
cts_raw_gene <- txi.sum_scaledTPM$counts          #  37 361 x 41  

# remove transcript/genes that have zeros across all the samples
cts_raw_transcript <- cts_raw_transcript[rowSums(cts_raw_transcript) > 0,]  # 151 271 x 41 
cts_raw_gene <- cts_raw_gene[rowSums(cts_raw_gene) > 0,]                    #  31 896 x 41

# save raw datasets
if(output_save==TRUE){
  write.csv(cts_raw_transcript[,-outliers_trans], "statistical_testing/input/counts_raw_transcript-level.csv")
  write.csv(cts_raw_gene[,-outliers_gene], "statistical_testing/input/counts_raw_gene-level.csv")
}

### => normalised <= ###
# NOTE: remove the outliers before normalisation (pre-normalised counts)
cts_pre_norm_transcript <- cts_raw_transcript[,-outliers_trans]
cts_pre_norm_gene <- cts_raw_gene[,-outliers_gene]
  
# perform normalisation using log-CPM (these counts can be used for visualisation purposes only)
# source: https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/cpm
# https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

cts_norm_transcript <- cpm(cts_pre_norm_transcript, log = TRUE)
cts_norm_gene <- cpm(cts_pre_norm_gene, log = TRUE)

# save raw datasets
if(output_save==TRUE){
  write.csv(cts_norm_transcript, "statistical_testing/input/counts_norm_transcript-level.csv")
  write.csv(cts_norm_gene, "statistical_testing/input/counts_norm_gene-level.csv")
}

# perform counts transformation using VST (these counts can be used for statisitcal testing)
# source: http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://rdrr.io/bioc/DESeq2/man/varianceStabilizingTransformation.html

# vst() requires integers as input
storage.mode(cts_pre_norm_transcript) <- "integer" 
storage.mode(cts_pre_norm_gene) <- "integer" 

vst_transcript <- vst(cts_pre_norm_transcript, blind=FALSE)
vst_gene <- vst(cts_pre_norm_gene, blind=FALSE)

# save VST transformed datasets
if(output_save==TRUE){
  write.csv(vst_transcript, "statistical_testing/input/counts_VST_transcript-level.csv")
  write.csv(vst_gene, "statistical_testing/input/counts_VST_gene-level.csv")
}

# other normalisation methods:
# there are many other normalisation methods, such as calcNormFactors function (edgeR package),
# source: http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# which require information about groups so it would need to be performed right before Mann-Whitney testing, for each feature
