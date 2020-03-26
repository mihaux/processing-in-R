# this script changes annotation from Entrez gene id to gene symbols  (output from featureCounts)


# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db"); library(org.Hs.eg.db)


# load dataset
df <- read.csv("/nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv", row.names=1)

df[1:3,1:3]

my_keys <- c("100287102", "653635", "102466751", "100302278")




#selecting
select(org.Hs.eg.db,
       keys = my_keys,
       columns=c("ENTREZID","SYMBOL","GENENAME"),
       keytype="ENTREZID")
