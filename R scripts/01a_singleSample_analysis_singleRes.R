##########################################################################
### Script to visualize cluster signals at a given cluster resolution
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
  res <- 0.7            # 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.8

# sub-directory creation #
  data.dir <- file.path(main.path,"data")
  sample.dir <- file.path(main.path,"results","single samples",sample.name)
  clus.dir<-file.path(sample.dir,paste0("clus@",res))
  if (!file.exists(clus.dir)){dir.create(clus.dir)}
  marker.file <- "markerList.txt"
  markers <- read.table(file.path(data.dir,marker.file),sep = '\t', header = T)
  
################################################################################
### Marker levels on slice and clusters at a given resolution
################################################################################
  dot.size <- 0.8
  codex.obj.clus<-readRDS(file.path(sample.dir,paste0(sample.name,".rds")))
  res_col<-paste0("Akoya_snn_res.",res) 
  # slide images and heatmap at a given resolution
  p1<- ImageDimPlot(codex.obj.clus, fov = sample.name, group.by = res_col, 
                   size=dot.size, border.color = "NA") +
    scale_y_reverse()+
    NoGrid()
  ggsave(file.path(clus.dir,paste0(sample.name,"_Image_@",res,".png")),
         p1, device = "png", width=12, height=10,)
  # cluster images and heatmap at a given resolution
  p1<- ImageDimPlot(codex.obj.clus, fov = sample.name, group.by = res_col, 
                    split.by = res_col, size=dot.size, border.color = "NA") +
    NoLegend()+
    scale_y_reverse()+
    NoGrid()
  p2<-DoHeatmap(codex.obj.clus, group.by = res_col,features = markers$MarkerName, size = 3) +
    NoLegend()
  ggsave(file.path(clus.dir,paste0(sample.name,"_ClusterImages_Heatmap_@",res,".png")),
         p1+p2, device = "png", width=20, height=10)
  # DotPlot of markers at given resolution
  pdf(file.path(clus.dir,paste0(sample.name,"_DotPlot_@",res,".pdf")),
      width=20, height=10,useDingbats=T)
  print(DotPlot(codex.obj.clus, features=markers$MarkerName,
                group.by=res_col,col.min = -0.4, col.max = 4,
                cols = c("white","red"),dot.scale=8))
  dev.off()
  # nCounts in clusters at given resolution
  pdf(file.path(clus.dir,paste0(sample.name,"_nCounts_@",res,".pdf")),
      width=20, height=10,useDingbats=T)
  print(VlnPlot(codex.obj.clus,features = "nCount_Akoya" ,group.by = res_col, pt.size = 0) +
          theme(plot.title = element_text(hjust = 0)))
  dev.off()
  # marker levels on slice image and clusters at a given resolution
  for (i in 1:length(markers$MarkerName)) {
    p1 <- ImageFeaturePlot(codex.obj.clus, fov = sample.name, 
                           features = markers$MarkerName[i],
                           min.cutoff = "q10", max.cutoff = "q90") +
      scale_y_reverse()
    p2 <- FeaturePlot(codex.obj.clus,features = markers$MarkerName[i],reduction = "umap") +
      xlab("UMAP 1") +
      ylab("UMAP 2")
    p3 <- DimPlot(codex.obj.clus, group.by=res_col, reduction="umap", label=T, label.size=6, pt.size=0.2) +
      ggtitle(paste("Clusters res",res))+
      xlab("UMAP 1") +
      ylab("UMAP 2") +
      theme(plot.title = element_text(hjust = 0)) +
      NoLegend()
    p4 <- VlnPlot(codex.obj.clus,features = markers$MarkerName[i],group.by = res_col, pt.size = 0, slot = "scale.data")
    ggsave(file.path(clus.dir,paste0(sample.name,"_",markers$MarkerName[i],"_@",res,".png")),
           p1+p2+p3+p4, device = "png", width=16, height=16)
  }
  
