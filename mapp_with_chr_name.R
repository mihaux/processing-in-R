
library(EnsDb.Hsapiens.v75)
library(stringr)

main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES_archived/run_15_Mar21/entire_output_directory"

temp = list.files(path=main_dir, pattern="quant.sf", recursive = TRUE, full.names = TRUE)

myfiles = lapply(temp, read.delim, stringsAsFactors=FALSE)

# all matrices has same transcript IDs, need to map them with chromosomal positions

edb <- EnsDb.Hsapiens.v75

ids <- str_sub(as.vector(myfiles[[1]]$Name), start = 1, end = 15)

## Get the data
mapped_chr <- select(edb, keys=ids, columns=c("SEQNAME"), keytype="TXID")

# exclude those transcripts that were not mapped
ids_new <- ids[ids %in% mapped_chr$TXID]

myfiles_new <- lapply(myfiles, function(x) x[ids %in% mapped_chr$TXID,])

# split in half: chrX+Y vs the rest
sub_XY <- ids_new[c(which(mapped_chr$SEQNAME == "X"), which(mapped_chr$SEQNAME == "Y"))]
sub_rest <- ids_new[-c(which(mapped_chr$SEQNAME == "X"), which(mapped_chr$SEQNAME == "Y"))]

# create new .sf file for 2 subsets:
myfiles_sub_XY <- lapply(myfiles_new, function(x) x[c(which(mapped_chr$SEQNAME == "X"), which(mapped_chr$SEQNAME == "Y")),])
myfiles_sub_rest <- lapply(myfiles_new, function(x) x[-c(which(mapped_chr$SEQNAME == "X"), which(mapped_chr$SEQNAME == "Y")),])

# get sample IDs
smp_names <- str_sub(temp, start = 109, end = 116)

dir_out_XY <- "/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES_archived/run_15_Mar21/output_sub_XY/"

for (i in seq(1:40)) {
  dir.create(paste0(dir_out_XY, smp_names[i]))
  filename = paste0(dir_out_XY, smp_names[i], "/", "quant.sf")
  write.table(myfiles_sub_XY[[i]], filename, row.names=FALSE, sep="\t")
}

dir_out_rest <- "/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES_archived/run_15_Mar21/output_sub_rest/"

for (i in seq(1:40)) {
  dir.create(paste0(dir_out_rest, smp_names[i]))
  filename = paste0(dir_out_rest, smp_names[i], "/", "quant.sf")
  write.table(myfiles_sub_rest[[i]], filename, row.names=FALSE, sep="\t")
}
