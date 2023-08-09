##########################################################################
### Script to annotate clusters with cell phenotypes
##########################################################################

library(Seurat)
library(ggplot2)

project <- "HNSCC"

### main settings ###
  sample.name <- "s2" #"s2","s3","s4","s6","s7","s10","s11","s12","s13"
  # do edit "folder1","folder2","..." with your folder path
  main.path <- file.path("folder1","folder2","...",project)
  res <- 0.7            # 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.8

# sub-directory creation #
  data.dir <- file.path(main.path,"data")
  sample.dir <- file.path(main.path,"results","single samples",sample.name)
  annot.dir<-file.path(sample.dir,paste0("clus_annot@",res))
  if (!file.exists(annot.dir)){dir.create(annot.dir)}
  marker.file <- "markerList.txt"

################################################################################
### Cluster re-labeling
################################################################################
  codex.obj.clus <- readRDS(file.path(sample.dir,paste0(sample.name,".rds")))
  clus.annot <- read.table(file.path(sample.dir,paste0(sample.name,"_clusterAnnot@",res,".txt")),
                           sep='\t',header = T,stringsAsFactors = F)
  # cluster re-labeling #
  codex.obj.clus[["single.annot"]]<-codex.obj.clus[[paste0("Akoya_snn_res.",res)]]
  codex.obj.clus[["merge.annot"]]<-codex.obj.clus[[paste0("Akoya_snn_res.",res)]]
  codex.single.annot<-codex.obj.clus[["single.annot"]]
  codex.merge.annot<-codex.obj.clus[["merge.annot"]]
  for (i in 1:dim(clus.annot)[1]) {
    levels(codex.single.annot$single.annot)[levels(codex.single.annot$single.annot) == clus.annot[i,1]] <- clus.annot[i,"annot"]
    levels(codex.merge.annot$merge.annot)[levels(codex.merge.annot$merge.annot) == clus.annot[i,1]] <- clus.annot[i,"merged"]
  }
  codex.obj.clus[["single.annot"]]<-codex.single.annot$single.annot
  codex.obj.clus[["merge.annot"]]<-codex.merge.annot$merge.annot

################################################################################
### Save Seurat object
################################################################################
  saveRDS(codex.obj.clus, file = file.path(sample.dir,paste0(sample.name,"_annoted.rds")))
  
################################################################################
### Images and plots with re-labelled clusters at given resolution
################################################################################
  dot.size <- 0.8
  markers <- read.table(file.path(data.dir,marker.file),sep = '\t', header = T)
  # dimensionality reductions with re-labeled clusters #
  pdf(file.path(annot.dir,paste0(sample.name,"_annot_tSNE_@",res,".pdf")),
      width=20, height=10,useDingbats=T)
  p1 <- DimPlot(codex.obj.clus, group.by = "single.annot",reduction = "tsne", label = T,
          label.size = 5, label.box = T) +
    ggtitle("Single cluster annotation")+
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    theme(plot.title = element_text(hjust = 0)) +
    NoLegend()
  p2 <- DimPlot(codex.obj.clus, group.by = "merge.annot",reduction = "tsne", label = T,
                label.size = 5, label.box = T) +
    ggtitle("Merged annotation")+
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    theme(plot.title = element_text(hjust = 0)) +
    NoLegend()
  print(p1+p2)
  dev.off()
  # slide images of any single cell population #
  pdf(file.path(annot.dir,paste0(sample.name,"_annot_clusters_@",res,".pdf")),
      width=12, height=10,useDingbats=T)
  print(ImageDimPlot(codex.obj.clus, fov = sample.name, group.by = "single.annot", 
                     split.by = "single.annot", size=dot.size, border.color = "NA") +
          scale_y_reverse()+
          NoGrid())
  print(ImageDimPlot(codex.obj.clus, fov = sample.name, group.by = "merge.annot", 
               split.by = "merge.annot", size=dot.size, border.color = "NA") +
          scale_y_reverse()+
          NoGrid())
  dev.off()
  # slide image with all cell populations #
  pdf(file.path(annot.dir,paste0(sample.name,"_annot_image_@",res,".pdf")),
      width=12, height=10,useDingbats=T)
  print(ImageDimPlot(codex.obj.clus, fov = sample.name, group.by = "single.annot", 
                     size=dot.size, border.color = "NA") +
          scale_y_reverse()+
          NoGrid())
  print(ImageDimPlot(codex.obj.clus, fov = sample.name, group.by = "merge.annot", 
                     size=dot.size, border.color = "NA") +
          scale_y_reverse()+
          NoGrid())
  dev.off()
  # marker heatmap in all cell populations #
  p1 <- DoHeatmap(codex.obj.clus, group.by = "merge.annot",features = markers$MarkerName, size = 3) +
    NoLegend()
  ggsave(file.path(annot.dir,paste0(sample.name,"_annot_heatmap_@",res,".png")),
         p1, device = "png", width=20, height=10)
 