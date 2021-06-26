# exploring Swish output files (using ID-11026 as example)

# define main directory
main_dir <- "/Users/ummz/Desktop/ID-11026"

# /aux_info
# cmd_info.json
# /libParams
# lib_format_counts.json
# /logs
# quant.sf

# load 'quant.sf' file 
quant_file <- read.delim(file.path(main_dir, "quant.sf"), stringsAsFactors = FALSE)

# drop suffixes
quant_file_names <- substr(quant_file$Name, 1, 15)

# load gene / transcript matrices
dir_gene <- "/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES_archived/swish_v2/ALL_SAMPLES_41/output_all/counts/gene_level_counts.csv"
dir_transcript <- "/Users/ummz/Documents/OneDrive - University of Leeds/ANALYSES_archived/swish_v2/ALL_SAMPLES_41/output_all/counts/transcript_level_tximeta_counts.csv"
cts_gene <- read.csv(dir_gene, row.names = 1, stringsAsFactors = FALSE, check.names = F) 
cts_transcript <- read.csv(dir_transcript, row.names = 1, stringsAsFactors = FALSE, check.names = F)

# extract rows for SSTR2 gene and its 2 transcripts

# transcripts:                         SSTR2-201              SSTR2-202 
# ID                              | ENST00000357585.4   | ENST00000579323.5   |
# transcript length in bp         | 7685                | 530                 |
# protein length in amino acids   | 369 aa              | No protein          |

quant_file[which(quant_file_names == "ENST00000357585"),]
#                       Name    Length    EffectiveLength      TPM          NumReads
# 135293    ENST00000357585.3   7683      7222.726             0.322086     28.609

quant_file[which(quant_file_names == "ENST00000579323"),]
#                       Name    Length    EffectiveLength      TPM          NumReads
# 135292    ENST00000579323.5   530       388                  0            0

cts_transcript[which(rownames(cts_transcript) == "ENST00000357585"),c(1:2)]
#                 ID-11026 ID-11028
# ENST00000357585   28.609   47.994

cts_transcript[which(rownames(cts_transcript) == "ENST00000579323"),c(1:2)]
#                 ID-11026 ID-11028
# ENST00000579323        0        0

# gene: ENSG00000180616 (ENSG00000180616.9 - the suffix can be different)
# gene length: 71,161,149-71,172,772 = 11623 bp

cts_gene[which(rownames(cts_gene) == "ENSG00000180616"),c(1:2)]
#                 ID-11026 ID-11028
# ENSG00000180616   28.609   47.994


# check another gene that have 2 transcripts and both are protein coding

# gene ISOC1 => ENSG00000066583
cts_gene[which(rownames(cts_gene) == "ENSG00000066583"),c(1:2)]
#                 ID-11026 ID-11028
# ENSG00000066583      107      152

# transcripts:                         ISOC1-201             ISOC1-202 
# ID                              | ENST00000173527.6   | ENST00000514194.5   |
# transcript length in bp         | 1942                | 581                 |
# protein length in amino acids   | 298aa               | 188aa	              |

quant_file[which(quant_file_names == "ENST00000173527"),]
#                    Name Length EffectiveLength      TPM NumReads
# 67876 ENST00000173527.5   1940        1483.832 3.120772   56.948

quant_file[which(quant_file_names == "ENST00000514194"),]
#                   Name Length EffectiveLength      TPM NumReads
# 67875 ENST00000514194.5    581         408.696 9.958256   50.052

cts_transcript[which(rownames(cts_transcript) == "ENST00000173527"),c(1:2)]
#                 ID-11026 ID-11028
# ENST00000173527   56.948  109.247

cts_transcript[which(rownames(cts_transcript) == "ENST00000514194"),c(1:2)]
#                 ID-11026 ID-11028
# ENST00000514194   50.052   42.753

# transcript1 + transcript2   = gene
# 56.948      + 50.052        = 107
# 109.247     + 42.753        = 152


#### test different normalisation methods ####
# source: https://github.com/crazyhottommy/RNA-seq-analysis/blob/master/salmon_kalliso_STAR_compare.md

countToTpm <- function( counts, effLen ) {
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

cts_transcript_TPM <- countToTpm(quant_file$NumReads, quant_file$EffectiveLength)

# TPM for gene level are just the sums of TPMs for all transcripts of a gene

library(EnsDb.Hsapiens.v86)
library(dplyr)

salmon_output = quant_file
tx2gene = transcripts(EnsDb.Hsapiens.v86, 
                      columns=c("tx_id", "gene_name", "gene_id"), 
                      return.type="DataFrame")

#gene_tpms <- salmon_output %>%
#  inner_join(tx2gene, by=c("Name"="tx_id")) %>%
#  group_by(gene_id, gene_name) %>%
#  summarise(TPM=sum(TPM))



countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}


# TODO: check if you can get gene length




