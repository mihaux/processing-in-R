# processing-in-R

## 1) FastQC_report.R

USAGE: `Rscript FastQC_report.R [directory-to-results] [directory-for-output]`

EXAMPLE:
```
Rscript FastQC_report.R /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_IV/1_quality_control/report /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_IV/1_quality_control/postprocessed
```

## 2) STAR_stats.R 

USAGE: `Rscript STAR_stats.R [directory-to-results] [directory-for-output]`

EXAMPLE:
```
Rscript STAR_stats.R /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_I_Nov19/4_alignement/bam /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_I_Nov19/4_alignement
```

## 3) featureCounts_merge.R

USAGE: `Rscript featureCounts_merge.R [directory-to-results] [directory-for-output]`

EXAMPLE:
```
Rscript featureCounts_merge.R /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_I /Users/ummz/OneDrive\ -\ University\ of\ Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_I
```

## 4) downstream.R

USAGE: `Rscript downstream.R [counts.csv] [annotation.csv] [directory-for-output]`

EXAMPLE: 
```
Rscript downstream.R /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/CiC_Clinical_data_FINAL.csv /Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_IV_Feb20/6_downstream_analysis/single-end
```

## 5) changeAnno.R

USAGE: `Rscript changeAnno.R [counts.csv] [directory-for-output]`

EXAMPLE:
for mode_I
```
Rscript changeAnno.R mode_I /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/single-end/processed/mode_I/counts_merged.csv /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/single-end/processed/mode_I/
```
and
```
Rscript changeAnno.R mode_I /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/paired-end/processed/mode_I/counts_merged.csv /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/paired-end/processed/mode_I/
```
for mode_II
```
Rscript changeAnno.R mode_II /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/all_counts_SE.csv /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/single-end/processed/mode_II/
```
and
```
Rscript changeAnno.R mode_II /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/paired-end/processed/mode_II/all_counts_PE.csv /nobackup/ummz/analyses/run_IV_Feb20/5_featureCounts/paired-end/processed/mode_II/
```

