# ScSpliceShiner
[![source](https://img.shields.io/badge/Source_code-support-blue.svg)](https://github.com//1027victor/ScSpliceShiner/tree/main/R)
## Introduction
## What is ScSpliceShiner
`SpliceShiner` is a shiny application based on the R package MARVEL, which aims to make single-cell alternative splicing analysis and gene expression analysis accessible to users through an interactive interface.

If you use `ScSpliceShiner`, please remember to cite its orignal paper:

Wen W X, Mead A J, Thongjuea S. MARVEL: an integrated alternative splicing analysis platform for single-cell RNA sequencing data[J]. Nucleic Acids Research, 2023, 51(5): e29-e29.


### Overview of splicing events

![](/inst/app/www/overview of splicing events.jpeg)

### Differential analysis

![](/inst/app/www/Differential analysis.jpe)

### Gene ontology analysis

![](/inst/app/www/Gene ontology analysis.jpeg)


### Dynamics analysis

![](/inst/app/www/Dynamics analysis.jpeg)


### Usage Demo
- Overview of splicing events
- Differential analysis
- Gene ontology analysis
- Dynamics analysis


## Installation
### Install ScSpliceShiner in Rstudio
Before installing this App, you will need to install some **dependent R packages** on your R.

```
install.packages("shiny")
install.packages("ComplexHeatmap")
install.packages("bslib")
install.packages("circlize")
install.packages("colourpicker") 
install.packages("DT")
install.packages("dplyr")
install.packages("data.table")
install.packages("shinyWidgets")
install.packages("shinycssloaders")
install.packages("shinyjs")
install.packages("highcharter")
install.packages("MARVEL")
install.packages("ggplot2")
install.packages("Matrix")
install.packages("plyr")
install.packages("scales")
install.packages("factoextra")
install.packages("FactoMineR")
install.packages("fitdistrplus")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("gtools")
install.packages("kSamples")
install.packages("pheatmap")
install.packages("reshape2")
install.packages("S4Vectors")
install.packages("scales")
install.packages("stringr")
install.packages("textclean")
install.packages("twosamples")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("AnnotationDbi")
BiocManager::install("Biostrings")
BiocManager::install("BSgenome")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
BiocManager::install("clusterProfiler")
BiocManager::install("GenomicRanges")
BiocManager::install("IRanges")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("phastCons100way.UCSC.hg38")

```

Once you have completed the installation of the dependencies, start downloading and installing the ScSpliceShiner.
```{r}
devtools::install_github("1027victor/ScSpliceShiner")
```
Once the installation is complete, run `ScSpliceShiner::run_app()` to open the APP page.






