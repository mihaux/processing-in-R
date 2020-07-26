# script to perform Mann-Whitney U test (also known as Wilcoxon rank sum test or Mann-Whitney test)
# 2 methods: wilcox.test() [built-in R function] and wilcoxTest() [GSALightning package]

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

# install (if necessary) and load package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GSALightning", quietly = TRUE)) BiocManager::install("GSALightning"); library(GSALightning)

# get working directory to recognise the machine
w_dir <- getwd()

# create a shortcut for the OneDrive directory where all files are stored
if(startsWith(w_dir, "/Users/michal")){           
  main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"    # on my mac
} else if (startsWith(w_dir, "/Users/ummz")) {    
  main_dir <- "/Users/ummz/OneDrive - University of Leeds"                # on uni mac    
} else {
  print("Unrecognised machine.")
}

#args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
       \n(1 - input) path to _x.csv file with count data, 
       \n(2 - annotation) path to _x.csv annotation file and
       \n(3 - output) path where output files should be stored", call.=FALSE)
}

# NOTE !!! : THERE MUST BE A "/" AT THE END OF ARGUMENT 3

args <- c(paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/DESeq2/norm_counts.csv"),
          paste0(main_dir, "/data/metadata/clinical_data/cic_clinical_data_v2_split/cic_clinical_data_v2_summary_ORDERED.csv"),
          paste0(main_dir, "/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/statistical_testing/"))

# Example of usage: 
# Rscript test_Mann_Whitney.R 

cat("Directories with data (IN): "); cat(args[1], sep="\n")
cat("Directory for results (OUT): "); cat(args[3], sep="\n")
setwd(args[3])

# load normalised counts from DESeq
df <- read.csv(args[1], row.names = 1)

# change data format to matrix and integer
dat <- as.matrix(df)   
storage.mode(dat) <- "integer"             # class: matrix | type: integer

# load annotation (clinical) data which will be used for coldata
anno <- read.csv(args[2], row.names = 1)

# add "ID_" to all rownames
rownames(anno) <- paste0("ID_", rownames(anno))

# define groups to be compared
#group = anno$gender..1.male..2.female.                     # group labels as gender
group = anno$visual.loss.at.BL..0.no..1.yes.                # group labels as visual loss

# NOTE: anno$visual.loss.ever...0.no..1.yes. is exactly the same

# RUN STATISTICAL TESTING

# HYPOTHISES: we want to know if the mean of group 1 differs from the mean of group 2.

# NOTE: no metter the comparison group_1 vs group_2 or group_2 vs group_1, 
# the results are the same, just the 'statistic' component differs (but all is proportional)

# 'alternative' argument => "two,sided" (default)
# if you want to test whether the mean of group_1 is less than the mean of group_2    => alternative = "less"
# if you want to test whether the mean of group_1 is greater than the mean of group_2 => alternative = "greater"

# INTERPRETATION:
# if the p-value of the test is less than the significance level alpha = 0.05. 
# we can conclude that the mean of group 1 is significantly different from the mean of group 2 with a p-value of [p-value =].

# initiate list for results
res <- list()

# iteration over all genes (as we perform statistical testing for each gene)
for(j in 1:nrow(dat)){
  
  # get data frame for one gene at a time
  temp = dat[j, ]
  
  # perform statistical testing CHANGE TO group==1 AND group==2 WHEN RUNNING FOR 'gender
  res[[j]] = wilcox.test(temp[group==0], 
                         temp[group==1], 
                         alternative = "two.sided",
                         exact = FALSE)               # to suppress the warning message saying that “cannot compute exact p-value with tie”
                                                      # it comes from the assumption of a Wilcoxon test that the responses are continuous. 
}
# might need to add as.numeric() if does not work

# the object 'res' contains the result of the test | length = 52239
#res[[1]]
# Wilcoxon rank sum test with continuity correction
# data:  temp[group == 2] and temp[group == 1]
# W = 254.5, p-value = 0.07999
# alternative hypothesis: true location shift is not equal to 0

# in the OUTPUT: (assume there are 2 groups: A and B)
# => W is a rank sum subtracted by a constant. 
# => if group A is the reference level in the factor variable group, the rank sum is computed for the data in group A
# => the value of W is the rank sum for group A subtracted by nA(nA+1)/2. 
# => the value of W is not important, as we perform the Wilcoxon Mann-Whitney test to determine if there is significant difference between the two groups. 
# => we are primarily interested in the p-value
# => by default wilcox.test() calculates an exact p-value if the samples contain less than 50 finite values and there are no ties. 
# => Otherwise, a normal approximation is used. R uses normal approximation to calculate the p-value.
  
# retrieve t-statistic for each gene,
all_stat = unlist(lapply(res, function(x) x$statistic))

# retrieve its corresponding p-value
all_pval = unlist(lapply(res, function(x) x$p.value))

# Multiplicity adjustment
# After obtaining the p-values, you can make multiplicity correction to the pvalue. 
all_adj.pval <- p.adjust(all_pval, "BH")    # Benjamini & Hochberg (1995) ("BH" or its alias "fdr")

# result table
# create a table as R data frame
result.table2 = data.frame(ID=rownames(dat), statistic=all_stat, pvalue=all_pval, fdr.pvalue=all_adj.pval)

# sort the table based on the order of the (adjusted) p-values
result.table2.sorted = result.table2[order(all_adj.pval),]

# and then show the top 10 (most) significant genes.
result.table2.sorted[1:10,]       # listing the top 10 genes

# select only significant genes
significant.cutoff <- length(which(result.table2.sorted$fdr.pvalue < 0.05))

result.table2.significant <- result.table2.sorted[1:significant.cutoff,]

# get the number of genes with pval < 5 % and padj < 5 %
length(which(result.table2$pvalue < 0.05))
length(which(result.table2$fdr.pvalue < 0.05))

# save results
#write.csv(result.table2.sorted, file="results_table_Mann-Whitney_gender_norm.csv")
write.csv(result.table2.sorted, file="results_table_Mann-Whitney_visual-loss_norm.csv")

#---------------------------- use another method GSALightning package ----------------------------#

# source: https://www.bioconductor.org/packages/release/bioc/vignettes/GSALightning/inst/doc/vignette.html

# the GSALightning package offers the Mann-Whitney U test for single-gene testing. 
# Mann-Whitney U test is the non-parametric version of the independent t-test for two-sample problem. 
# To perform the Mann-Whitney U test, call the wilcoxTest() function:
  
#singleWilcox <- wilcoxTest(eset = dat, fac = factor(anno$gender..1.male..2.female.), tests = "unpaired")
singleWilcox <- wilcoxTest(eset = dat, fac = factor(anno$visual.loss.at.BL..0.no..1.yes.), tests = "unpaired")


# save results
#write.csv(singleWilcox, file="results_table_Mann-Whitney_gender_GSALightning_norm.csv")
write.csv(singleWilcox, file="results_table_Mann-Whitney_visual-loss_GSALightning_norm.csv")








