# this script changes annotation from Entrez gene id to gene symbols  (output from featureCounts)

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler"); library(clusterProfiler)

args <- commandArgs(trailingOnly = TRUE)

# files to be processed:

# /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/single-end/processed/mode_I/counts_merged.csv
# /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/paired-end/processed/mode_I/counts_merged.csv

# /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv
# /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/paired-end/processed/mode_II/all_counts_PE.csv

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
	\n(1 - mode) either 'mode_I' or 'mode_II'
	\n(2 - input) path to .csv file with count data
	\nand (3 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[2], sep="\n")

cat("Directory for results (OUT): ")
cat(args[3], sep="\n")

setwd(args[3])

# load count data

if(args[1] == "mode_I"){

	df <- read.csv(args[2], row.names = 2)
	
	# create object containing rownames
	my_key <- rownames(df)
	head(my_key)

	# mapping between keys and columns
	df_out <- select(org.Hs.eg.db, keys = my_key, columns=c("ENTREZID","SYMBOL","GENENAME"), keytype="ENTREZID")
 
	# write output table
	write.csv(df_out, file="mapped_EntrezID_mode_I.csv")

} else if (args[1] == "mode_II"){

	df <- read.csv(args[2], row.names = 1)
	
	# create object containing rownames
	my_key <- rownames(df)
	head(my_key)

	# mapping between keys and columns
	df_out <- select(org.Hs.eg.db, keys = my_key, columns=c("ENTREZID","SYMBOL","GENENAME"), keytype="ENTREZID")

	# write output table
	write.csv(df_out, file="mapped_EntrezID_mode_II.csv")

} else {
	stop("ERROR: running mode parameter must be defined as either mode_I or mode_II", call.= FALSE)
}

############################################################
########## using Bitr - Biological Id TRanslator #########
# Usage bitr(geneID, fromType, toType, OrgDb, drop = TRUE)
############################################################
# Ian's data
dir_2 <- "/Users/ummz/Documents/OneDrive - University of Leeds/comparison_with_Ian_results/featureCounts_IAN/"

# read count tables
counts_2_dups <- read.csv(paste(dir_2, "gene_counts_withDups.csv", sep = ""), row.names = 1)
counts_2_nodups <- read.csv(paste(dir_2, "gene_counts_withoutDups.csv", sep = ""), row.names = 1)








