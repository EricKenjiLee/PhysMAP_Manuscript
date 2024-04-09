# PhysMAP

R code for multi-modal analysis to delineate cell types using electrophysiology.

The following dependencies are not a strict requirement but are those in which we know are compatible. This project is <ins>not</ins> compatible with Seurat v5 as there were breaking changes introduced into our usage of anchor alignment.

:warning:**WARNING**:warning:: If you are a reviewer of this manuscript, ___please do not star or fork this repo!___ Doing so breaks reviewer anonymity.

Note: While this code reproduces figures in the associated manuscript, system (operating system) differences will result in plots that differ in appearance but this does not materially affect any results (see manuscript Supp. Fig. 2).

## Generating Manuscript Figures

Regenerating the figures/results found in the manuscript Lee _et al._ 2024 can be done by simply running the scripts described below for each session. The only requirement is to first install Seurat v4 (and specifically v4) as newer versions of Seurat are incompatible with this analysis. To do this, the following steps should be followed.

### Installation of Seurat v4 (Hao and Hao _et al._ 2021)

### Mouse S1 Juxtacellular Dataset (Yu _et al._ 2019)

### Mouse A1 Probe Extracellular Dataset (Lakunina _et al._ 2020)

### Mouse Visual Cortex and Hippocampus Datasets (Petersen _et al._ 2021)

### Cell Type Classifier (Petersen _et al._ 2021)

## Session Info
R version 4.2.2 (2022-10-31)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Ventura 13.5.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
aricode_1.0.3      R.matlab_3.7.0     cowplot_1.1.3      SeuratObject_5.0.1 Seurat_4.1.0       here_1.0.1         reshape2_1.4.4     nnet_7.3-18       caret_6.0-94       lattice_0.20-45    lubridate_1.9.3    forcats_1.0.0     stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2        readr_2.1.5       tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.0      tidyverse_2.0.0  
