##########################################################################
### Script to generate the Seurat object from raw .CSV files
### and normalize and cluster the data in the UMAP space 
### Each .CSV file has a "cys:" before any marker,
### has "cell_id:cell_id" in the first column, and has
### duplicated labels for any other column with a ":" as "region:region"
##########################################################################

library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(clustree)
library(openxlsx)

project <- "HNSCC"

### main settings ###
  sample.name <- "s2" #"s2","s3","s4","s6","s7","s10","s11","s12","s13"
  # do edit "folder1","folder2","..." with your folder path
  main.path <- file.path("folder1","folder2","...",project)

# sub-directory creation #
  data.dir <- file.path(main.path,"data")
  upper.result.dir<-file.path(main.path,"results")
  if (!file.exists(upper.result.dir)){dir.create(upper.result.dir)}
  main.result.dir<-file.path(upper.result.dir,"single samples")
  if (!file.exists(main.result.dir)){dir.create(main.result.dir)}
  result.dir <- file.path(main.result.dir,sample.name)
  if (!file.exists(result.dir)){dir.create(result.dir)}
  marker.file <- "markerList.txt"
  markers <- read.table(file.path(data.dir,marker.file),sep = '\t', header = T)
  data.filename <- paste0(sample.name,".csv")
  
####################################
### Data pre-processing 
####################################
  # load Akoya CODEX data #
  codex.obj <- LoadAkoya(filename = file.path(data.dir,data.filename),
                         type = "processor",fov = sample.name)
  # normalization and scaling #
  codex.obj <- NormalizeData(object = codex.obj, normalization.method = "CLR", margin = 2)
  codex.obj <- ScaleData(codex.obj)
  pdf(file.path(result.dir,paste0("01_",sample.name,"_Norm_total_expression_after+before_norm.pdf"))) 
  par(mfrow = c(3,1))
  hist(colSums(as.matrix(codex.obj@assays$Akoya@counts)), breaks=100, 
       main=paste0(sample.name, " - Total expression before normalization"), xlab="Sum of expression") 
  hist(colSums(as.matrix(codex.obj@assays$Akoya@data)), breaks=100, 
       main=paste0(sample.name, " - Total expression after normalization"), xlab="Sum of expression")
  hist(colSums(as.matrix(codex.obj@assays$Akoya@scale.data)), breaks=100, 
       main=paste0(sample.name, " - Total expression after normalization and scaling"), xlab="Sum of expression")
  dev.off()
  
################################################################################
### Determine statistically significant principal components 
################################################################################
  n.PCs <- 20
  VariableFeatures(codex.obj) <- rownames(codex.obj)  # since the panel is small, treat all features as variable.
  codex.obj.temp <- RunPCA(object = codex.obj, npcs = n.PCs, verbose = T)
  # 1. exploring PCs to determine relevant sources of heterogeneity 
  pdf(file.path(result.dir,paste0("02_",sample.name,"_Heatmap.pdf")))
  DimHeatmap(codex.obj.temp, dims=1:9, cells = 1000, balanced = TRUE)
  DimHeatmap(codex.obj.temp, dims=10:20, cells = 1000, balanced = TRUE)
  dev.off()
  # 2. plot of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph  
  pdf(file.path(result.dir,paste0("03_",sample.name,"_PCElbowPlot.pdf")))
  ElbowPlot(object=codex.obj.temp, ndims=n.PCs)
  dev.off()
  
###########################################################
### Dimensionality reduction (PCA, UMAP, t-SNE)
###########################################################
  # Do edit the number of PCs to be used
  pca_dims <- 11
  codex.obj <- RunPCA(object = codex.obj, npcs = pca_dims, verbose = TRUE)
  codex.obj <- RunTSNE(object = codex.obj, dims = 1:pca_dims, verbose = TRUE)
  codex.obj <- RunUMAP(object = codex.obj, dims = 1:pca_dims, verbose = TRUE)
  point.size=0.2
  dim.reduct <- c("PCA","tSNE","UMAP")
  plot <- NULL
  for (i in 1:length(dim.reduct)) {
    if (dim.reduct[i]=="PCA"){prefix <- "PC"}else{prefix <- dim.reduct[i]}
    plot[[i]] <- DimPlot(codex.obj, reduction=tolower(dim.reduct[i]), pt.size=point.size) +
      ggtitle(paste0(sample.name," - ",dim.reduct[i])) +
      theme(plot.title = element_text(hjust = 0))+
      theme(legend.position="none")+
      xlab(paste(prefix,"1")) +
      ylab(paste(prefix,"2"))
  }
  pdf(file.path(result.dir, paste0("04_", sample.name, "_Embeddings.pdf")), width=2*7, height=2*7, useDingbats=FALSE)
  print(plot_grid(plotlist=plot, ncol=2, nrow=2))
  invisible(dev.off())

###########################################################
### clustering
###########################################################
  res.try <- seq(0.2,0.8,by=0.1)
  codex.obj <- FindNeighbors(object = codex.obj, dims = 1:pca_dims, verbose = T)
  codex.obj <- FindClusters(object = codex.obj, resolution = res.try)
  # cluster dynamics along resolutions  
  pdf(file.path(result.dir,paste0("05a_",sample.name,"_clusters_tree.pdf")), width=14, height=10)
  metadata<-codex.obj@meta.data
  for (res in res.try) {
    res_col<-paste0("CODEX_snn_res.",res)
    metadata[,res_col]<-as.integer(as.character(metadata[,res_col]))
  }
  clustree(metadata[,grep("CODEX_snn",names(metadata))], prefix = "CODEX_snn_res.")
  dev.off()
  # clusters on UMAP embedding
  n.col <- 3
  n.row <- ceiling(length(res.try)/n.col)
  plot <- NULL
  for (r in res.try){
    i <- match(r,res.try)
    res_col <- paste0("CODEX_snn_res.",r)
    plot[[i]] <- DimPlot(codex.obj, group.by=res_col, reduction="umap", 
                         label=T, label.size=6, pt.size=point.size) +
    ggtitle(paste("Clusters res",r))+
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(plot.title = element_text(hjust = 0)) +
    NoLegend()
  }
  pdf(file.path(result.dir, paste0("06a_", sample.name, "_clusters_UMAP.pdf")), width=n.col*7, height=n.row*7, useDingbats=FALSE)
  print(plot_grid(plotlist=plot, ncol=n.col, nrow=n.row))
  invisible(dev.off())
  # clusters on t-SNE embedding
  n.col <- 3
  n.row <- ceiling(length(res.try)/n.col)
  plot <- NULL
  for (r in res.try){
    i <- match(r,res.try)
    res_col <- paste0("CODEX_snn_res.",r)
    plot[[i]] <- DimPlot(codex.obj, group.by=res_col, reduction="tsne", 
                         label=T, label.size=6, pt.size=point.size) +
      ggtitle(paste("Clusters res",r))+
      xlab("t-SNE 1") +
      ylab("t-SNE 2") +
      theme(plot.title = element_text(hjust = 0)) +
      NoLegend()
  }
  pdf(file.path(result.dir, paste0("06b_", sample.name, "_clusters_tSNE.pdf")), width=n.col*7, height=n.row*7, useDingbats=FALSE)
  print(plot_grid(plotlist=plot, ncol=n.col, nrow=n.row))
  invisible(dev.off())
  
################################################################################
### Image and marker plots for clusters at all resolutions
################################################################################
  dot.size <- 0.8
  # Clusters on slide image at all resolutions
  pdf(file.path(result.dir,paste0("07_",sample.name,"_clusters_Images.pdf")),useDingbats=T,width=12, height=8)
  for (res in res.try) {
    res_col<-paste0("CODEX_snn_res.",res)
    print(ImageDimPlot(codex.obj, fov = sample.name, size=dot.size, group.by = res_col, border.color = "NA") +
            NoGrid() +
            scale_y_reverse())
    print(ImageDimPlot(codex.obj, fov = sample.name, group.by = res_col, 
                       split.by = res_col, size=dot.size, border.color = "NA") +
            NoLegend()+
            scale_y_reverse()+
            NoGrid())
  }
  dev.off()
  # Clusters on slide image and heatmap at all resolutions
  pdf(file.path(result.dir,paste0("08_",sample.name,"_clusters_Images_Heatmaps.pdf")),width=20, height=10, useDingbats=T)
  for (res in res.try) {
    res_col<-paste0("CODEX_snn_res.",res)
    p1 <- ImageDimPlot(codex.obj, fov = sample.name, group.by = res_col,
                       split.by = res_col, size=dot.size, border.color = "NA") +
      NoLegend() +
      scale_y_reverse() +
      NoGrid()
    p2 <-DoHeatmap(codex.obj, group.by = res_col,features = markers$MarkerName, size = 3) +
      NoLegend()
    print(p1+p2)
  }
  dev.off()
  # DotPlot of markers on clusters at all resolutions
  pdf(file.path(result.dir,paste0("09_",sample.name,"_clusters_DotPlot.pdf")),
      width=20, height=10,useDingbats=T)
  for (res in res.try) {
    res_col<-paste0("CODEX_snn_res.",res)
    print(DotPlot(codex.obj, features=markers$MarkerName,
                  group.by=res_col,col.min = -0.4, col.max = 4,
                  cols = c("white","red"),dot.scale=8) +
            ylab(paste0("cluster @ res ",res)) +
            NoLegend())
  }
  dev.off()

################################################################################
### Export average marker levels in clusters at all resolutions
################################################################################
  # Average marker levels in clusters at all resolutions
  wb<-createWorkbook(title="results")
  bold.style <- createStyle(textDecoration = "Bold")
  for (res in res.try) {
    res_col<-paste0("CODEX_snn_res.",res)
    clus.scale.ave <- AverageExpression(codex.obj,
                                        features = markers$MarkerName,
                                        group.by = paste0("CODEX_snn_res.",res),
                                        slot = "scale.data")
    clus.ave <- AverageExpression(codex.obj,
                                  features = markers$MarkerName,
                                  group.by = paste0("CODEX_snn_res.",res),
                                  slot = "data")
    addWorksheet(wb,paste0("res ",res))
    modifyBaseFont(wb, fontSize = 14, fontName = "Arial")
    colnames(clus.scale.ave$Akoya) <- paste0("c.",colnames(clus.scale.ave$Akoya))
    colnames(clus.ave$Akoya) <- paste0("c.",colnames(clus.ave$Akoya))
    writeData(wb, sheet = paste0("res ",res), round(clus.scale.ave$Akoya,3),
              rowNames = T, headerStyle = bold.style)
    conditionalFormatting(wb, paste0("res ",res),
                          cols = 1:ncol(clus.scale.ave$Akoya)+1, rows = 1:nrow(clus.scale.ave$Akoya)+1,
                          style = c("steelblue3", "white","indianred2"),
                          type = "colourScale")
    start.Col <- ncol(clus.scale.ave$Akoya)+3
    end.Col <- ncol(clus.scale.ave$Akoya)+3+ncol(clus.ave$Akoya)+1
    writeData(wb, sheet = paste0("res ",res), round(clus.ave$Akoya,3), startCol = start.Col,
              rowNames = T, headerStyle = bold.style)
    conditionalFormatting(wb, paste0("res ",res),
                          cols = start.Col:end.Col, rows = 1:nrow(clus.ave$Akoya)+1,
                          style = c("steelblue3", "white","indianred2"),
                          type = "colourScale")
  }
  saveWorkbook(wb, file.path(result.dir,paste0("10_",sample.name,"_clusters_averages.xlsx")), overwrite = TRUE)    

################################################################################
### Save Seurat object
################################################################################
  saveRDS(codex.obj, file = file.path(result.dir,paste0(sample.name,".rds")))
  
