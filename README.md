# BF591-Final_Project

#### Data Visualization (GSE64810)

### Data

Data used to build this R Shiny app was extracted from Gene Expression Omnibus (GEO): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810

This website uses the Huntington's Disease dataset from (Labadorf et al., 2016) titled 'mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington's Disease and neurologically normal individuals'")

### Tab Information

1. The first tab on this website contains information about the Sample metadata which can be accessed from `data/sample_metadata.csv` under this repository
2. The second tab uses the counts matrix to understand the count structure and aid in gene filtering for visualising Diagnostic plots, Heatmap and PCA. The data can be accessed from `data/norm_counts.csv` under this repository.
3. The third tab uses results from Differential Expression Analysis and generates a summary table alongwith volcano plots. The data can be accessed from `data/deseq_diff_exp_res.csv` under this repository.
4. The fourth tab uses Normalized Counts and Sample metadata for Individual Gene Expression Visualization. The data can be accessed from `data/sample_metadata.csv` and `data/norm_counts.csv` under this repository.

#### Reference: `Labadorf, A., Hoss, A., Lagomarsino, V., Latourelle, J., Hadzi, T., Bregu, J., MacDonald, M., Gusella, J., Chen, J., Akbarian, S., Weng, Z. and Myers, R., 2016. Correction: RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression. PLOS ONE, 11(7), p.e0160295.`
