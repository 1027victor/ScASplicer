# ScASplicer
[![source](https://img.shields.io/badge/Source_code-support-blue.svg)](https://github.com//1027victor/ScSpliceShiner/tree/main/R)
![Version](https://img.shields.io/badge/version-1.1.0-blue.svg)
![License](https://img.shields.io/github/license/1027victor/ScSpliceShiner.svg)
![Languages](https://img.shields.io/github/languages/top/1027victor/ScSpliceShiner.svg)
![Stars](https://img.shields.io/github/stars/1027victor/ScSpliceShiner.svg)
![Forks](https://img.shields.io/github/forks/1027victor/ScSpliceShiner.svg)
![Contributors](https://img.shields.io/github/contributors/1027victor/ScSpliceShiner.svg)

## Introduction
## What is ScASplicer
`ScASplicer` is a  R package that extends MARVEL's capabilities to support multiple cell populations, offering an interactive, code-free platform for alternative splicing and gene expression analysis. These advancements significantly improve MARVEL’s usability and functionality, enabling broader and more complex single-cell AS studies( main support Single-cell-RNA-sequencing data generated from plate-based library preparation methods such as Smart-seq2).

![](/inst/app/www/pipeline.jpg)

If you use `ScASplicer`, please remember to cite its orignal paper:

Wen W X, Mead A J, Thongjuea S. MARVEL: an integrated alternative splicing analysis platform for single-cell RNA sequencing data[J]. Nucleic Acids Research, 2023, 51(5): e29-e29.

## Features

### Overview of splicing events

![](/inst/app/www/overview_of_splicing_events.jpeg)

### Modality analysis

![](/inst/app/www/Modality_analysis.jpeg)

### Differential analysis

![](/inst/app/www/Differential_analysis.jpeg)

### Gene ontology analysis

![](/inst/app/www/Gene_ontology_analysis.jpeg)


### Dynamics analysis

![](/inst/app/www/Dynamics_analysis.jpeg)

## Installation
1. Install the R [(LINK)](https://cran.r-project.org/)
2. Install the free version of rStudio [(LINK)](https://www.rstudio.com/products/rstudio/download/)
3. Run the following command in rStudio to install ScASplicer as an R package:
```{r,eval=FALSE}
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
install.packages("MARVEL")
devtools::install_github("1027victor/ScASplicer")
```
Once the installation is complete, run the following command in rStudio to open the APP page.
```
ScASplicer::run_app()
```
