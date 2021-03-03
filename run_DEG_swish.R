






# load counts (raw)
cts_gene <- read.csv(paste0(data_dir, "counts_gene-level_no_outliers.csv"), row.names = 1)
cts_transcript <- read.csv(paste0(data_dir, "counts_transcript-level_no_outliers.csv"), row.names = 1)

# define path to annotation files (clinical / histological)
dir_anno <- paste0(main_dir, "/data/metadata/outliers_excluded/cic_clinical_data_v2_summary_ORDERED_outliers_excluded.csv")
#dir_anno <- paste0(main_dir, "/data/metadata/outliers_excluded/slide_scores_v6_outliers_excluded.csv")

# load annotation (clinical) data which will be used for coldata
anno <- read.csv(dir_anno, row.names = 1)


# create a list-based data object to store the counts
y_gene <- SummarizedExperiment(assays=cts_gene, colData=anno)                      # dim = 37361 x 40
y_transcript <- SummarizedExperiment(assays=cts_transcript, colData=colnames(cts_transcript))    # dim = 175775 x 40

###### Differential transcript expression ######

# (1) Running Swish at the transcript level

# Running swish has three steps: 
# => scaling the inferential replicates, (NOT NEEDED, USED DATA GENERATED WITH length-scaledTPM)
# => labeling the rows with sufficient counts for running differential expression, 
# => and then calculating the statistics. 

# NOTE: a random seed needs to be set before running swish(), as swish makes use of pseudo-random number generation in breaking ties and in calculating permutations

# The default number of permutations in swish is nperms=100. (leave it as it is)

library(fishpond)

#y <- scaleInfReps(y_gene)
y <- labelKeep(y_gene)

#y <- scaleInfReps(y)
y <- labelKeep(y)
y <- y[mcols(y)$keep,]
set.seed(1)
y <- swish(y_gene, x="gender", pair="line", nperms=64)



#### EXAMPLE ####
makeInfReps# make sim data for 30 samples:
y <- makeSimSwishData(n=30)

y$batch <- factor(rep(rep(1:3,each=5),2))
y <- scaleInfReps(y)
y <- labelKeep(y)
y <- swish(y, x="condition", cov="batch")


yy <- makeInfReps(y_gene, numReps, minDisp = 0.001)
library(fishpond)
