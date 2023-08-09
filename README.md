# CODEX_HNSCC
R scripts accompanying our manuscript "The tumor immune microenvironment architecture correlates with risk of recurrence in head and neck squamous cell carcinoma".

# Table of Contents

- [Introduction](https://github.com/bicciatolab/CODEX_HNSCC#introduction)
- [System requirements](https://github.com/bicciatolab/CODEX_HNSCC#system-requirements)
- [Scripts](https://github.com/bicciatolab/CODEX_HNSCC#scripts)

## Introduction
This repository comprises the scripts used to analyze the CODEX data and generate results presented in our publication “The tumor immune microenvironment architecture correlates with risk of recurrence in head and neck squamous cell carcinoma".  

All scripts are contained in the “R scripts” directory and are implemented in R based on Seurat functions for the [analysis of CODEX spatial data](https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#human-lymph-node-akoya-codex-system). Some experience with installing R packages will be required to get up and running. Only a familiarity with R is required to follow and customize the scripts for own use.

Key features of these scripts include:

1. Creation of the Seurat data object
2. Data preprocessing and cell phenotyping
2. Neighborhood identification
3. Quantification of pairwise cell-cell contacts

Single sample raw data and the single-cell data table of clustered, annotated cell types with metadata can be downloaded from [Mendeley](https://data.mendeley.com/datasets/t2yvtwnjx7).

#### Contact:

silvio.bicciato@unimore.it

#### Citation:

Weed DT, Zilio S, McGee C, Marnissi B, Sargi Z, Franzmann E, Thomas G, Leibowitz J, Nicolli E, Arnold D, Bicciato S, Serafini P. The tumor immune microenvironment architecture correlates with risk of recurrence in head and neck squamous cell carcinoma, _Cancer Research_ (2023)

## System requirements
* R version: >= 4.0.0
* Dependencies: *circlize*, *ClusterR*, *clustree*, *corrplot*, *cowplot*, *data.table*, *dplyr*, *ggplot2*, *ggrepel*, *gplots*, *grid*, *gridExtra*, *openxlsx*, *paletteer*, *pheatmap*, *RANN*, *RColorBrewer*, *REAT*, *RTriangle*, *scales*, and *Seurat*.

All scripts expect a project directory (e.g., "HNSCC") with a nested "data" directory containing all single sample aw data as .CSV files:

```
  Project directory
   |_ data

```
This directory structure can be easily customized in the "main settings" and "sub-directory" sections of each script.


## Scripts
1.	`01_singleSample_analysis.R`: this scripts preprocesses the raw .CSV data files using Seurat functions to [analyze image-based spatial data](https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#human-lymph-node-akoya-codex-system) in R >=4.0.0. Briefly, the CSV data file of each single sample from CODEX MAV is loaded in R using Seurat `LoadAkoya` function, and protein signals are normalized with the centered log-ratio based normalization. To detect the cell phenotypes, `01_singleSample_analysis.R` first performs Principal Component Analysis (PCA) for dimensionality reduction and then cluster analysis in the low-dimensional space. Before applying PCA, the data are scaled and the number of components for downstream analysis is determined. Then, graph-based clustering at different resolutions (e.g., ranging from 0.2 to 0.8 in steps of 0.1) is performed. Cell clusters are visualized on a protein intensity-based uniform manifold approximation and projection (UMAP) and on their spatial location. To associate cell phenotypes to clusters, levels of protein markers are displayed as dot plots for each cluster at any explored clustering resolution.
2. `01a_singleSample_analysis_singleRes.R`: for each sample, the optimal number of clusters (i.e., the clustering resolution) is determined based on visual inspection of cluster location and marker expression. Clusters with a similar morphological appearance in the tissue and similar marker expression profiles can be merged, and artifacts removed.
3. `02_singleSample_CellTypeAnnotation.R`: this script assigns cell type annotations to clusters at the selected resolution based on the average expression of protein markers in each cluster and visual inspection. The script uses a text file `sample.name_clusterAnnot@0.8.txt` listing the single and merged cell phenotype assigned to each cluster. A `sample.name_clusterAnnot@res.txt` file has to be created for each sample and placed in the "sample.name" sub-directory inside the "results" directory:

 ```
  Project directory
   |_ results
     |_sample.name

 ```
 
4. `03_Neighborhood_Identification.R`: this script identifies cellular neighborhoods (CNs), i.e., regions with a characteristic local composition of cell phenotypes, as presented in _Schürch CM, Bhate SS, Barlow GL, Phillips DJ, Noti L, Zlobec I, et al. Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front. Cell 2020;182:1341-59.e19_. For each cell of any sample, the script identifies a window consisting of W nearest neighboring cells (including the center cell) using the `nn2` function of the `RANN` R package (version 2.6.1). The `nn2` function uses a k-dimensional tree to find a given number of near neighbors (here, W) for each point identified by the X and Y coordinates in the input dataset. These windows are then clustered by their composition with respect to the cell types previously identified by graph-based clustering and supervised annotation. Specifically, each window is converted to a vector containing the frequency of each cell type among the W neighbors. Subsequently, windows are clustered using the `MiniBatchKmeans` function of the `ClusterR` package (version 1.2.9) implementing the Mini-batch K-means clustering algorithm with a given value of K. Each cell is then allocated to the CN of its surrounding window using the `predict_MBatchKMeans` function of the `ClusterR` package.
5. `04_CellCellContacts_w/inCN_single.R`: this script quantifies the pairwise cell-cell contacts within CNs. It starts from determining the direct neighbors of each cell within each cellular neighborhood using the Delaunay triangulation implemented in the `triangulate` function of the `RTriangle` R package (version 1.6-0.11). Then, it calculates the number of contacts between cells of types i and j as the set of edges Nijk between cells of type i and j in cellular neighborhood k returned by the triangulation. Finally, it displays pairwise contacts among different phenotype cells as circle plots using the `chordDiagram` function of the `circlize` R package (version 0.4.15).
6. `04a_CellCellContacts_w/inCN_metrics.R`: this script is used to evaluate if cell-cell contacts of cells assigned to a given phenotype are influenced by the cellular neighborhood. It compares the relative cell-cell contact frequencies (RF) for all phenotypes in each pair of CNs (e.g., TLS1 and TLS2) across all patients containing both CNs. The script calculates the relative cell-cell contact frequency between cells of types i and j in each cellular neighborhood k as Nijk/Nik where Nijk is the number of edges of the triangulation between cells of type I and j and Nik is the number of edges of type i in CN k. Finally, it compares the relative cell-cell contact frequency in pairs of CNs using the `ttest` function of the R `stats` package after discrding cell-cell interactions between pairs of phenotypes present in less than three samples.
