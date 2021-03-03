# script to perform Mann-Whitney U test (also known as Wilcoxon rank sum test or Mann-Whitney test)
# using the wilcox.test() built-in R function
# for HISTOLOGICAL variables
# this file was created based on "test_Mann-Whitney.R", which includes more detais about the test

#----------------------------------------------------------------------------------------------------#
# very clear source: http://courses.atlas.illinois.edu/spring2016/STAT/STAT200/RProgramming/NonParametricStats.html

# The unpaired two-samples Wilcoxon test (also known as Wilcoxon rank sum test or Mann-Whitney test) 
# => is a non-parametric alternative to the unpaired two-samples t-test, 
# => can be used to compare two independent groups of samples. 
# => is used when your data are not normally distributed.

# NOTE: if there are confounders, then you would need to rather use something like logistic regression, not Mann-Whitney test

# WHAT TO TEST: statistical test between 2 groups: visual loss vs no visual loss,
# to see if patients with visual loss are more likely to have GCA or not (if not, then there's no association)

# To perform two-samples Wilcoxon test comparing the means of two independent samples (x & y), 
# we perform a t-test on each gene (i.e. each row) by running the function wilcox.test() for each row, in a two-sample setting between group_1 and group_2.
#----------------------------------------------------------------------------------------------------#

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

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

args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=5) {
  stop("5 arguments must be supplied: 
       \n(1 - input) path to _x.csv file with count data (or just the directory), 
       \n(2 - annotation) path to _x.csv annotation file,
       \n(3 - feature) name of the feature for running and
       \n(4 - run_ID) running ID (i.e. type of data to be used)
       \n(5 - output) path where output files should be stored", call.=FALSE)
}

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 5

if(FALSE){
args <- c(paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/statistical_testing/input/counts_VST_transcript-level.csv"),
          paste0(main_dir, "/data/metadata/outliers_excluded/slide_scores_v6_outliers_excluded.csv"), 
          "GCA_present", 
          "VST_transcript-level",
          paste0(main_dir, "/ANALYSES_archived/run_13_Jan21/statistical_testing/output/histological/"))
}
# args[4] = {VST_transcript-level, VST_gene-level, raw_transcript-level, raw_gene-level}

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[5], sep="\n")
setwd(args[5])

#####################################################################################################
# load counts matrix
dds <- read.csv(args[1], row.names = 1, header = TRUE)

# change data format to matrix (required for wilcox.test() function later on)
dat <- as.matrix(dds)

# define running ID
run_id <- args[4]

# load annotation (clinical) data which will be used for coldata
anno <- read.csv(args[2], row.names = 1)

# NOTE: sample "ID_14058" is missing in 'slide_sc', needs to be removed from counts matrix as well
# if running for slide scores => remove sample "ID.14058" from all three datasets

dat <- dat[,-which(colnames(dat) == "ID.14058")]

# modify Occlusion.grade. and Intima.pattern. to make it 2 comparison groups
Occlusion_grade_new <- c(1:39)
Occlusion_grade_new[which(anno$Occlusion_grade < 3)] <- 0      # less than 3           => replace with 0
Occlusion_grade_new[which(anno$Occlusion_grade >= 3)] <- 1     # equal or more than 3  => replace with 1
anno$Occlusion_grade <- Occlusion_grade_new

Intima_pattern_new <- c(1:39)
Intima_pattern_new[which(anno$Intima_pattern <= 1)] <- 0       # 0 and 1                => replace with 0
Intima_pattern_new[which(anno$Intima_pattern > 1)] <- 1        # 2 and 3                => replace with 1
anno$Intima_pattern <- Intima_pattern_new
  
Media_pattern_new <- c(1:39)
Media_pattern_new[which(anno$Media_pattern <= 1)] <- 0       # 0 and 1                => replace with 0
Media_pattern_new[which(anno$Media_pattern > 1)] <- 1        # 2 and 3                => replace with 1
anno$Media_pattern <- Media_pattern_new

# get index of the feature for running
running <- args[3]
run_ind <- which(colnames(anno) == running)

# define groups for comparison
group = as.vector(unlist(anno[run_ind]))

# initiate a list for results
res <- list()

# iteration over all genes (as we perform statistical testing for each gene)
for(j in 1:nrow(dat)){
  temp = dat[j, ]
  
  # perform statistical testing
  res[[j]] = wilcox.test(temp[group==0], temp[group==1], alternative = "two.sided", exact = FALSE)               
            # to suppress the warning message saying that “cannot compute exact p-value with tie”
            # it comes from the assumption of a Wilcoxon test that the responses are continuous. 
}

# retrieve p-value for each transcript/gene
all_pval = unlist(lapply(res, function(x) x$p.value))

# Multiplicity adjustment; Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
all_adj.pval <- p.adjust(all_pval, "fdr") 

# create a table as R data frame
result.table2 = data.frame(ID=rownames(dat), pvalue=all_pval, fdr.pvalue=all_adj.pval)

# sort the table based on the order of the p-values
result.table2.sorted_pval = result.table2[order(all_pval),]

# sort the table based on the order of the adjusted p-values
result.table2.sorted_padj = result.table2[order(all_adj.pval),]

# save tables
write.csv2(result.table2.sorted_pval, file=paste0("table_sorted_by_pvalues_", running, "_", run_id, ".csv"))
write.csv2(result.table2.sorted_padj, file=paste0("table_sorted_by_padjusted_", running, "_", run_id, ".csv"))
