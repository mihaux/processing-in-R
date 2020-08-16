# PCA tutorial; source: ???

# Principal component analysis (PCA) is one of the fundamental tools of population genetics for identifying sample clustering and outliers
# As explained in the talk by Gavin PCA is way of reducing a multidimensional space of data points to a set of orthogonal principal components that represent "directions" of greatest variance.  
# These may reflect population structure, but might also reflect cryptic relationships, poor QC, or other effects that lead to differences in genotype distribution.

############################## in PLINK ############################## 
# https://www.cog-genomics.org/
  
#We will use plink
# ->	to remove closely-related samples
# ->	to compute principal components
# ->	and to compute the SNP weights or loadings that tell us how principal components are weighted across the genome.

# create new directory and set it as working directory
dir.create("/Users/ummz/PCA_tutorial")
setwd("/Users/ummz/PCA_tutorial")

# A note on quality control

# Before carrying out a genetic analysis, like PCA, it's important to have a good-quality dataset, and this typically means carrying out careful quality control (QC) first.  

# The cleaned data consists of genotype calls at different sites.

### LD pruning of SNPs

# Genetic drift and other processes lead to linkage disequilibrium (LD) between SNPs along a chromosome.
# To ensure the PCs we compute represent genome-wide structure (not local LD) we'll first carry out LD pruning of our SNP set.
# This removes correlated pairs of SNPs so that the remaining SNPs are roughly independent.  (It also makes subsequent computations quicker.)  

# Run the following command to prune the dataset:
### plink --vcf chr19-clean.vcf.gz --maf 0.01 --indep-pairwise 50 5 0.2 --out chr19-clean 

# The above command tells plink to load the file chr19-clean.vcf.gz and to prune SNPs to leave SNPs with MAF at least 1%, with no pairs remaining with r2>0.2.  
# (The other parameters, here 50 and 5, affect how the computation works in windows across the genome.  
# You can read about the behaviour here: http://www.cog-genomics.org/plink2/ld).

Q. Look at the screen output from the above plink command.  How many variants were in the original dataset?  How many were removed because their frequency was below 1%?  How many variants were removed due to LD pruning?  How many variants remain?
  
  Type ls or use the file manager to view the directory.  The command above produced a number of files that all begin with the chr19-clean prefix.  For our purposes, the most important one is chr19-clean.prune.in, as this lists the SNPs that remain after pruning.  Feel free to look at these files using less or a text editor.






