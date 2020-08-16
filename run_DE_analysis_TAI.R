
#---------------------------------------------------------------------------------------------------
# this package does not seem to by completly suitable for my analyses
#---------------------------------------------------------------------------------------------------

# use myTAI package: https://cran.r-project.org/web/packages/myTAI/vignettes/Expression.html
#install.packages("myTAI", build_vignettes = TRUE, dependencies = TRUE); library(myTAI)

# Some of the methods to to detect differentially expressed genes are based on non-statistical quantification of expression differences (e.g. fold-change and log-fold-change), 
# but most methods are based on statistical tests to quantify the significance of differences in gene expression between samples. 
# These statistical methods can furthermore be divided into two methodological categories: parametric tests and non-parametric tests. 
# The DiffGenes() function available in myTAI implements the most popular and useful methods to detect differentially expressed genes.

# NOTE: when using DiffGenes() it is assumed that your input dataset has been normalized before passing it to DiffGenes(). 
# For RNA-Seq data DiffGenes() assumes that the libraries have been normalized to have the same size, 
# i.e., to have the same expected column sum under the null hypothesis (or the lib.size argument in DiffGenes() is specified accordingly).

# Fold-Changes

# A fold change in gene expression is simply the ratio of the gene expression level of one sample against a second sample: 
# ‘Ei1 / Ei2’ , where 
# ‘Ei1’ is the expression level of gene ‘i’ in sample one and 
# ‘Ei2’ is the expression level of gene 'i’ in sample two.

# EXAMPLE:
#data("PhyloExpressionSetExample")

# Detection of DEGs using the fold-change measure
#DEGs <- DiffGenes(ExpressionSet = PhyloExpressionSetExample[1:5,1:8],
#                  nrep          = 2,
#                  method        = "foldchange",
#                  stage.names   = c("S1","S2","S3"))

#head(DEGs)

# When selecting method = "log-foldchange" it is assumed that the input ExpressionSet stores log2 expression levels. 
# transform absolute expression levels to log2 expression levels using the tf() function before log-fold-changes are computed.


# Wilcoxon-Mann-Whitney test (Mann-Whitney U test) -> non-parametric (Welch t-test is parametric)

# The Wilcoxon-Mann-Whitney test is a nonparametric test to quantify the shift in empirical distribution parameters. 
# Nonparametric tests are useful when sample populations do not meet the test assumptions of parametric tests.

# Assumptions about input data
# => independent samples
# => continuous data
# => (approximate) normality

# Nevertheless, although in most cases log2 expression levels are used to perform the Welch t-test assuming that expression levels are log-normal distributed which approximates a normal distribution in infinity, in most cases the small number of replicates is not sufficient enough to fulfill the (approximate) normality assumption of the Welch t-test.
# Due to this fact, non-parametric, sampling based, or generalized linear model based methods have been proposed to quantify p-values of differential expression. Nevertheless, the DiffGenes() function implements the Welch t-test for the detection of differentially expressed genes, allowing users to compare the results with more recent DEG detection methods/methodologies also implemented in DiffGenes().

# load normalised counts from DESeq
#counts_norm <- read.csv("/Users/michal/Documents/OneDrive - University of Leeds/ANALYSES/downstream/rerun_FINAL_July20/rerun_5/DESeq2/norm_counts.csv", row.names = 1)

# Detection of DEGs using the p-value returned by a Wilcoxon-Mann-Whitney test
#Wilcox.DEGs <- DiffGenes(ExpressionSet = counts_norm,
#                         method        = "wilcox.test",
#                         stage.names   = c("S1","S2","S3"))

# add p.adjust.method = "BH" for p-value correction

# look at the results
#Wilcox.DEGs

# adjust p-values by specifying the p.adjust.method argument.

# Detection of DEGs using the p-value returned by a Wilcoxon-Mann-Whitney test
# and furthermore, adjust p-values for multiple comparison
# using the Benjamini & Hochberg (1995) method: method = "BH"
# and filter for significantly differentially expressed genes (alpha = 0.05)