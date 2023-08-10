# CODEX_HNSCC
R scripts accompanying our manuscript "The tumor immune microenvironment architecture correlates with risk of recurrence in head and neck squamous cell carcinoma".

# Table of Contents

- [Introduction](https://github.com/bicciatolab/CODEX_HNSCC#introduction)
- [System requirements](https://github.com/bicciatolab/CODEX_HNSCC#system-requirements)
- [Scripts](https://github.com/bicciatolab/CODEX_HNSCC#scripts)

## Introduction
This repository contains scripts used to analyze CODEX data and generate the results presented in our publication titled "The architecture of the tumor immune microenvironment correlates with the risk of recurrence in head and neck squamous cell carcinoma." 

All scripts are located within the "R scripts" directory and are implemented in R, utilizing Seurat functions for [analysis of CODEX spatial data](https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#human-lymph-node-akoya-codex-system). Some experience with installing R packages will be necessary to set up and use these scripts. However, only a basic familiarity with R is needed to understand and customize the scripts for personal use.

Key features of these scripts comprise:

1. Generation of the Seurat data object
2. Preprocessing of data and identification of cell phenotypes
3. Identification of cellular neighborhoods
4. Quantification of pairwise cell-cell contacts

The raw data for individual samples and the single-cell data table containing clustered and annotated cell types, along with associated metadata, can be obtained from [Mendeley](https://data.mendeley.com/datasets/t2yvtwnjx7).

#### Contact:

silvio.bicciato@unimore.it, pserafini@miami.edu

#### Citation:

Weed DT, Zilio S, McGee C, Marnissi B, Sargi Z, Franzmann E, Thomas G, Leibowitz J, Nicolli E, Arnold D, Bicciato S, Serafini P. The tumor immune microenvironment architecture correlates with risk of recurrence in head and neck squamous cell carcinoma, _Cancer Research_ (2023)

## System requirements
* R version: >= 4.0.0
* Dependencies: *circlize*, *ClusterR*, *clustree*, *corrplot*, *cowplot*, *data.table*, *dplyr*, *ggplot2*, *ggrepel*, *gplots*, *grid*, *gridExtra*, *openxlsx*, *paletteer*, *pheatmap*, *RANN*, *RColorBrewer*, *REAT*, *RTriangle*, *scales*, and *Seurat*.

All scripts require a project directory (e.g., "HNSCC") with an embedded "data" directory containing all single sample in the form of .CSV files:

```
  Project directory
   |_ data
```
This directory structure can be easily tailored within the "main settings" and "sub-directory" segments of each script.


## Scripts
1.	`01_singleSample_analysis.R`: this script preprocesses the raw .CSV data files utilizing Seurat functions to [analyze image-based spatial data](https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#human-lymph-node-akoya-codex-system) in R version 4.0.0 or later. Briefly, each CSV data file representing a single sample from CODEX MAV is loaded into R using the Seurat `LoadAkoya` function, followed by normalization of protein signals through centered log-ratio based normalization.
For cell phenotype identification, the `01_singleSample_analysis.R` script initially performs Principal Component Analysis (PCA) for dimensionality reduction, subsequently conducts cluster analysis within the reduced-dimensional space. Prior to PCA, the data are scaled and the number of components for downstream analysis is determined. Then, graph-based clustering is performed at multiple resolutions (e.g., ranging from 0.2 to 0.8 in increments of 0.1). Cell clusters are visualized on a protein intensity-based uniform manifold approximation and projection (UMAP), as well as on their spatial location. To associate cell phenotypes with clusters, levels of protein markers are displayed as dot plots for each cluster at any explored clustering resolution.
2. `01a_singleSample_analysis_singleRes.R`: for each sample, the optimal number of clusters (referred to as clustering resolution) is determined through visual assessment of cluster distribution and marker expression. Clusters exhibiting comparable morphological characteristics within the tissue and similar marker expression profiles can be considered for merging, and artifacts can be eliminated.
3. `02_singleSample_CellTypeAnnotation.R`: this script assigns cell type annotations to clusters at the selected resolution based on the average expression of protein markers within each cluster and visual examination. The script uses a text file named `sample.name_clusterAnnot@0.8.txt` which lists the individual and merged cell phenotypes assigned to each cluster. For each sample, a `sample.name_clusterAnnot@res.txt` file must be created and placed in the "sample.name" sub-directory inside the "results" directory:

 ```
  Project directory
   |_ results
    |_sample.name
 ```
 
4. `03_Neighborhood_Identification.R`: this script detects cellular neighborhoods (CNs), referring to regions characterized by a distinctive local composition of cell phenotypes, as demonstrated in the work by _Schürch CM, Bhate SS, Barlow GL, Phillips DJ, Noti L, Zlobec I, et al. Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front. Cell 2020;182:1341-59.e19_. For each cell within any sample, the script identifies a window comprising the W nearest neighboring cells, inclusive of the central cell. This operation is executed by the `nn2` function from the R package `RANN` (version 2.6.1). By employing a k-dimensional tree, the `nn2` function identifies the specified number of close neighbors (W) for every point defined by the X and Y coordinates in the input dataset. Subsequently, these windows are grouped based on their composition in relation to the cell types earlier identified through graph-based clustering and supervised annotation. More specifically, each window is transformed into a vector containing the occurrence frequency of each cell type among the W neighbors. Following this, the windows are subjected to clustering using the `MiniBatchKmeans` function from the `ClusterR` package (version 1.2.9), which implements the Mini-batch K-means clustering algorithm with a designated value of K. Finally, each cell is allocated to the CN of its surrounding window employing the `predict_MBatchKMeans` function from the `ClusterR` package.
5. `04_CellCellContacts_w/inCN_single.R`: this script quantifies the interactions between pairs of cells within CNs. The process begins by identifying the immediate neighbors of each cell within each cellular neighborhood. This is accomplished through the application of Delaunay triangulation, which is executed using the `triangulate` function from the R package `RTriangle` (version 1.6-0.11). Then, the script computes the count of contacts between cells of types i and j. This calculation is carried out based on the set of edges Nijk connecting cells of type i and j within cellular neighborhood k, as determined by the triangulation process. In the final step, the script displays these pairwise contacts among cells of various phenotypes in the form of circle plots using the `chordDiagram` function from the R package `circlize` (version 0.4.15).
6. `04a_CellCellContacts_w/inCN_metrics.R`: this script is used to assess whether the cell-cell interactions among cells of a specific phenotype are influenced by their cellular neighborhood. It conducts a comparison of the relative frequencies (RF) of cell-cell contacts  for all phenotypes within each pair of CNs (such as TLS1 and TLS2), across all patients containing both CNs. In the process, the script computes the relative cell-cell contact frequency between cells of types i and j within each cellular neighborhood k. This calculation is accomplished by dividing Nijk (the count of edges in the triangulation connecting cells of type i and j) by Nik (the count of edges of type i in CN k). Finally, it compares the relative cell-cell contact frequencies across pairs of CNs. This is carried out using the `ttest` function from the R package `stats` on cell-cell interactions between pairs of phenotypes present in at least three samples.
