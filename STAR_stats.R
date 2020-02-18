# STAR processing script that allows to merge all [sample_ID]_Log.final.out files and generate a summary table

# TODO:

# install (if necessary) and load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# args <- commandArgs(trailingOnly = TRUE)

args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_I_Nov19/4_alignement/bam", 
          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_I_Nov19/4_alignement")

if (length(args)!=2) {
  stop("2 arguments must be supplied: \n(1 - input) path to directory with data and \n(2 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[1], sep="\n")

cat("Directory for results (OUT): ")
cat(args[2], sep="\n")

setwd(args[1])
temp = list.files(pattern="*.out")
all_files = lapply(temp, function(x) read.csv(x,sep = "\t", header = FALSE))

setwd(args[2])

# extract IDs from file names 
IDs <- sub("_.*", "", temp)

# add colnames to each table
for (i in 1:length(all_files)) {
  names(all_files[[i]]) <- c("V1", IDs[i])
}

# merge all all_files
merged <- do.call("cbind", all_files)

# remove every second column starting from the 3rd one
to_delete <- seq(3, length(merged), 2)
merged_final <- merged[,-to_delete]

#Started job on                                   | Dec 12 18:46:27
#Started mapping on                               | Dec 12 18:46:57
#Finished on                                      | Dec 12 18:55:53
#Mapping speed, Million of reads per hour         |          144.29
#Number of input reads                            |        21483273
#Average input read length                        |             150
#UNIQUE READS:                
#Uniquely mapped reads number                     |        17943309
#Uniquely mapped reads %                          |          83.52%
#Average mapped length                            |          150.44
#Number of splices: Total                         |         9798337
#Number of splices: Annotated (sjdb)              |         9709065
#Number of splices: GT/AG                         |         9714551
#Number of splices: GC/AG                         |           59166
#Number of splices: AT/AC                         |            6910
#Number of splices: Non-canonical                 |           17710
#Mismatch rate per base, %                        |           0.51%
#Deletion rate per base                           |           0.04%
#Deletion average length                          |            1.90
#Insertion rate per base                          |           0.00%
#Insertion average length                         |            1.31
#MULTI-MAPPING READS:                
#Number of reads mapped to multiple loci          |         2282825
#% of reads mapped to multiple loci               |          10.63%
#Number of reads mapped to too many loci          |          444554
#% of reads mapped to too many loci               |           2.07%
#UNMAPPED READS:                
#Number of reads unmapped: too many mismatches    |               0
#% of reads unmapped: too many mismatches         |           0.00%
#Number of reads unmapped: too short              |          633819
#% of reads unmapped: too short                   |           2.95%
#Number of reads unmapped: other                  |          178766
#% of reads unmapped: other                       |           0.83%
#CHIMERIC READS:                
#Number of chimeric reads                         |               0
#% of chimeric reads                              |           0.00%

# generate a table 


