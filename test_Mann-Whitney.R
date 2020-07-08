# script to perform Mann-Whitney U test

# http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r

# there was a link sent by Mark with some sources about this test
# other sources about statistical testing:
# https://www.mayo.edu/research/documents/parametric-and-nonparametric-demystifying-the-terms/doc-20408960
# https://www.ajodo.org/action/showPdf?pii=S0889-5406%2815%2900840-9

# The comparison groups would be: affected (case) vs unaffected (controls) for the GCA.present feature 
# and similarly for the other outcomes: visual loss vs. no visual loss for the visual.loss.at.BL feature, etc.

# the Mann–Whitney test will be used to see if there's a correlation between visual loss and presents of GCA; 
# look at genes one-by-one to see ; 
# the comparison groups are: visual loss vs no visual loss

# exposure: visual loss
# outcome: GCA presence

# WHAT TO TEST: statistical test between 2 groups: visual loss vs no visual loss,
# to see if patients with visual loss are more likely to have GCA or not (if not, then there's no association)

# CAN DO SIMALARLLY FOR OTHER CLINICAL FEATURES AS WELL




