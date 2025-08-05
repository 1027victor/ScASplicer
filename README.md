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

## 📂 Example Data and Tutorial Videos

To help users get started quickly, we provide one **example dataset** and two **instructional videos**. The example dataset can be used to test and reproduce the analysis workflow, while the videos demonstrate how to operate the tool step by step.

### ▶️ Tutorial Videos

| No. | Title                       | Link                                                                 |
|-----|-----------------------------|----------------------------------------------------------------------|
| 1   | MARVEL_pipeline             | [Watch Video](https://drive.google.com/file/d/1gw2FhMuj3E-bgaM4u7vAHVz_rpdxSwPO/view)       |
| 2   | Complete_Example_Demo.mp4   | [Watch Video] https://drive.usercontent.google.com/download?id=1SId4RsuEmDj_rMtzERodPOIECL5eksD7         |

> 💡 **Tip**: You can open the `.mp4` files directly in your browser or download and play them using any media player.

### 📁 test data

📥 **Download or explore the example dataset**: [test_data](https://drive.usercontent.google.com/download?id=1vdKI2qk54rOTuPBPNv-eVFterGy4P2DL)




## Features

### Overview of splicing events

![](/inst/app/www/overview_of_splicing_events.jpg)

### Modality analysis

![](/inst/app/www/Modality_analysis.jpg)

### Differential analysis

![](/inst/app/www/Differential_analysis.jpg)

### Gene ontology analysis

![](/inst/app/www/Gene_ontology_analysis.jpg)


### Dynamic analysis

![](/inst/app/www/Dynamics_analysis.jpg)

### Visualization of annotation tracks
![](/inst/app/www/visualization_of_annotation_tracks.jpeg)

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
BiocManager::install("trackViewer")
BiocManager::install("rtracklayer")
BiocManager::install("VariantAnnotation")
BiocManager::install("Rsamtools")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
install.packages("MARVEL")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("1027victor/ScASplicer")
```
Once the installation is complete, run the following command in rStudio to open the APP page.
```
ScASplicer::run_app()
```
## Install packages from conda
+ create conda env 
```
git clone https://github.com/1027victor/ScASplicer.git
cd ScASplicer
conda env create -f ScASplicer.yaml
```
+ activate env
```
conda activate ScASplicer
```
+ Run the following command in the R console
```
install.packages("MARVEL")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
devtools::install_github("1027victor/ScASplicer")
```
+ Then run the following command in the conda environment terminal
```
Rscript -e "options(browser = function(...) NULL, shiny.port = 10027, shiny.host = '0.0.0.0',sass.cache =FALSE); ScASplicer::run_app()"
```

Note:
+ If you use windows system ,The application will be run locally at http://127.0.0.1:10027, user could open the address with google chrome or other modern browsers.
+ If you use  Linux server,The application will be run locally at http://xxx.xxxx.xxx.xxx:10027, xxx.xxx.xxx.xxx is the IP address of the server, user could open the address with google chrome or other modern browsers.
## Use docker image

We have prepared docker images for `ScASplicer`. With docker installed, user could simplely invoke
the app with command below, and will be able to invoke the application directly from 
**"Containers/ Apps"** menu when opening Docker Desktop next time.

Pull the pre-built image from [dockerhub](https://hub.docker.com/), use:
```
docker pull biovictor520zy/scasplicer:victor
```
Run ScASplicer Docker Container with Port Mapping
```
docker run -itd --name scasplicer_hpw -p 10027:10027 biovictor520zy/scasplicer:victor
```
The application will be run locally at `http://xxxx.xxxx.xxxx.xxxx:10027`, user could open
the address with browsers.

Note:
+ If you use windows system ,The application will be run locally at http://127.0.0.1:10027, user could open the address with google chrome or other modern browsers.
+ If you use  Linux server,The application will be run locally at http://xxx.xxxx.xxx.xxx:10027, xxx.xxx.xxx.xxx is the IP address of the server, user could open the address with google chrome or other modern browsers.
- Windows users needs to install docker desktop, and type the same command
above in any terminal app, e.g. `PowerShell`.
- The image is a bit large, please reserve 5 GB space for it.  

