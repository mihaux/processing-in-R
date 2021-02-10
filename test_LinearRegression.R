
# You can test for correlation between the histological features (e.g. using Spearman's) but also you've shown signifcant association between STR1 and various histological features, but these are all univariate tests. 

# Instead regress STR1 (using linear regression) on each in turn, then assuming GCA_present is the most sig regress STR1 on GCA_present + each of the other hist features in turn (STR1 ~ GCA_present + Giant_cells, STR1 ~ GCA_present +intima_pattern etc.). 

# Then you can check the signifcance of each of the other features (i.e. is the coefficient sig different to 0?). It may be that there's no association with the other features beyond just GCA_present, so when this is included in the model, none of the other hist features are significant. 

# So you'll see whether the association between STR1 and each of these is effectively multiple interesting associations, or really just driven by one underlying association.


