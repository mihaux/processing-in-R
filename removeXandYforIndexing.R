# remove chrX and chrY from the fasta transcriptome file that is used to create an index for Salmon alignment

##-------------- in terminal --------------##

## unzip the compressed file
# gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

## extract every other line, starting from 1st
# awk 'NR%2==1' test.txt > test_extracted.txt

## extract the third column which contains the information about chromosome location
# cut -d " " -f 3 Homo_sapiens.GRCh38.cdna.all.fa > third_col.txt


##----------------------------------------##
setwd("/Users/ummz/Desktop/for_indexing")

# load the whole file
res <- read.delim("Homo_sapiens.GRCh38.cdna.all.fa", header = FALSE, stringsAsFactors = FALSE)
#tescik <- read.delim("test.txt", header = FALSE, stringsAsFactors = FALSE)

# add an extra column 
res$V2 <- seq(1:nrow(res))

# select the rows that starts with ">" 
to_extract <- which(startsWith(res$V1, ">"))
id_seq <- which(!startsWith(res$V1, ">"))

res_extracted <- res[to_extract,]

#------------------------------------------------------------#
# get the position to be removed in the original hyge dataset
res_split_all_1 <- strsplit(res$V1, split = "gene:ENSG")
res_split_all_2 <- unlist(lapply(res_split_all_1, function(x) x[1]))
res_split_all_3 <- strsplit(res_split_all_2, split = "cdna")
res_split_all_4 <- unlist(lapply(res_split_all_3, function(x) x[2]))
res_split_all_5 <- gsub(" ", "", res_split_all_4, fixed = TRUE)
res_split_all_6 <- strsplit(res_split_all_5, ":")
res_split_all_7 <- unlist(lapply(res_split_all_6, function(x) x[3]))

# for chrX => 5638 lines to be removed
length(which(res_split_all_7 == "X"))                   # 5609
length(which(res_split_all_7 == "CHR_HSCHRX_2_CTG3"))   # 12
length(which(res_split_all_7 == "CHR_HSCHRX_2_CTG12"))  # 5
length(which(res_split_all_7 == "CHR_HSCHRX_1_CTG3"))   # 12
length(which(res_split_all_7 == "Y"))                   # 630

unique(res_split_all_7)     # 380 (381 including one NA)

#------------------------------------------------------------#
res_split_1 <- strsplit(res_extracted$V1, split = "gene:ENSG")
#res_split_2 <- unlist(lapply(res_split_1, function(x) x[1]))
#res_split_3 <- strsplit(res_split_2, split = "cdna")
#res_split_4 <- unlist(lapply(res_split_3, function(x) x[2]))
#res_split_5 <- gsub(" ", "", res_split_4, fixed = TRUE)
#res_split_6 <- strsplit(res_split_5, ":")
#res_split_7 <- unlist(lapply(res_split_6, function(x) x[3]))

#sink("unique.txt")
#unique(res_split_7)
#sink()

# for chrX => 5638 lines to be removed
#which(res_split_7 == "X")                   # 5609
#which(res_split_7 == "CHR_HSCHRX_2_CTG3")   # 12
#which(res_split_7 == "CHR_HSCHRX_2_CTG12")  # 5
#which(res_split_7 == "CHR_HSCHRX_1_CTG3")   # 12
#which(res_split_7 == "Y")                   # 630

###########################################################################
# source: https://www.r-bloggers.com/2020/01/comparing-ensembl-gtf-and-cdna/
###########################################################################

# Download cDNA fasta file
library(here)
library(Biostrings)
library(stringr)
library(dplyr)
library(tidyr)
library(biomartr)

# Download cDNA fasta file
if (!file.exists(here("hs_cdna94.fa.gz"))) {
  download.file("ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
                destfile = here("hs_cdna94.fa.gz"))
}

cdna <- readDNAStringSet(here("hs_cdna94.fa.gz"))

# Extract transcript ID from fasta sequence name
cdna_tx <- str_extract(names(cdna), "^ENST\\d+")

# Extract annotation from fasta sequence names
cdna_meta <- tibble(transcript_id = cdna_tx,
                    cr = str_extract(names(cdna),"(?<=((chromosome)|(scaffold)):GRCh38:).*?(?=\\s)"),
                    gene_biotype = str_extract(names(cdna), "(?<=gene_biotype:).*?(?=\\s)"),
                    gene_id = str_extract(names(cdna), "(?<=gene:).*?(?=\\.)"),
                    gene_symbol = str_extract(names(cdna), "(?<=gene_symbol:).*?(?=\\s)")) %>% 
                        separate(cr, into = c("seqnames", "start", "end", "strand"), sep = ":") %>% 
                            mutate(start = as.integer(start), 
                                   end = as.integer(end), 
                                   strand = case_when(strand == "1" ~ "+", strand == "-1" ~ "-", TRUE ~ "*"))

# create subsets:
# "X"                   # 5609
# "CHR_HSCHRX_2_CTG3"   # 12
# "CHR_HSCHRX_2_CTG12"  # 5
# "CHR_HSCHRX_1_CTG3"   # 12
# "Y"                   # 630

length(which(cdna_meta$seqnames == "Y"))
length(which(cdna_meta$seqnames == "X"))
length(which(cdna_meta$seqnames == "CHR_HSCHRX_2_CTG3"))
length(which(cdna_meta$seqnames == "CHR_HSCHRX_2_CTG12"))
length(which(cdna_meta$seqnames == "CHR_HSCHRX_1_CTG3"))


only_y <- cdna[which(cdna_meta$seqnames == "Y")]

only_x <- cdna[c(which(cdna_meta$seqnames == "X"),
                 which(cdna_meta$seqnames == "CHR_HSCHRX_2_CTG3"),
                 which(cdna_meta$seqnames == "CHR_HSCHRX_2_CTG12"),
                 which(cdna_meta$seqnames == "CHR_HSCHRX_1_CTG3"))]

only_x_y <- cdna[c(which(cdna_meta$seqnames == "Y"),
                   which(cdna_meta$seqnames == "X"),
                   which(cdna_meta$seqnames == "CHR_HSCHRX_2_CTG3"),
                   which(cdna_meta$seqnames == "CHR_HSCHRX_2_CTG12"),
                   which(cdna_meta$seqnames == "CHR_HSCHRX_1_CTG3"))]

only_1_22 <- cdna[-c(which(cdna_meta$seqnames == "Y"),
                    which(cdna_meta$seqnames == "X"),
                    which(cdna_meta$seqnames == "CHR_HSCHRX_2_CTG3"),
                    which(cdna_meta$seqnames == "CHR_HSCHRX_2_CTG12"),
                    which(cdna_meta$seqnames == "CHR_HSCHRX_1_CTG3"))]

# write subset of FASTA transcript
writeFasta(only_y, "out_subsets/only_chrY.fa", mode = "a")
writeFasta(only_x, "out_subsets/only_chrX.fa", mode = "a")
writeFasta(only_x_y, "out_subsets/only_chrX-Y.fa", mode = "a")
writeFasta(only_1_22, "out_subsets/only_chr1-22.fa", mode = "a")

