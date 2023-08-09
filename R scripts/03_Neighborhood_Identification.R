########################################################################################
### Script to identify cell neighborhoods through neighborhood analysis
########################################################################################

library(Seurat)
library(RANN)
library(ClusterR)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggrepel)
library(corrplot)
library(openxlsx)
library(cowplot)
 
project <- "HNSCC"

### main settings ###
  sample.names <- c("s2","s3","s4","s6","s7","s10","s11","s12","s13")
  # do edit "folder1","folder2","..." with your folder path
  main.path <- file.path("folder1","folder2","...",project)

### CN settings ###
  cell.neigh <- 14   # number of total cell neighbors CN c(5,7,10)
  cell.window <- 10  # number of neighboring cells c(5,7,10)

### distance settings ###
  check.dist <- FALSE
  pixel.to.micron <- 0.377474 # conversion factor from pixel to micron
  min.dist <- 50 # minimum distance between cells in micron 

# sub-directory creation #
  data.dir <- file.path(main.path,"results","merged samples")
  main.result.dir<-file.path(main.path,"results","neighborhood analysis")
  if (!file.exists(main.result.dir)){dir.create(main.result.dir)}
  result.dir <- file.path(main.result.dir,paste0("CN",cell.neigh,"_window",cell.window))
  if (!file.exists(result.dir)){dir.create(result.dir)}
  markers <- read.table(paste0(file.path(main.path,"data"),"/",marker.file),sep='\t',header = T)
  marker.file <- "markerList.txt"
  
################################################################################
### Load Seurat objects for all samples
################################################################################
  obj.merged <- readRDS(file.path(data.dir,"merged.rds"))

######################################################################################
### Identification of cell.neigh CNs composed of cell.window neighboring cells
######################################################################################
  # construct the file for neighborhood analysis
    metadata <- data.frame(rownames(obj.merged[[]]),obj.merged[[]])
    colnames(metadata)[1] <- "EventID"
    keep.cols <- c("EventID","orig.ident","cell_id","region","tile_num","x","y","x_tile","y_tile", "size","merge.annot")
    dummy.matrix <- matrix(data = 0, nrow =dim(metadata)[1], ncol = length(unique(metadata$merge.annot)))
    colnames(dummy.matrix) <- unique(metadata$merge.annot)
    tmp.CN.input <- cbind(metadata[,keep.cols],dummy.matrix)
  # compute annotation of the cell.window neighboring cells in each single sample
    CN.input <- NULL
    for (i in 1 : length(sample.names)){
      sub.tmp.CN.input <- tmp.CN.input[tmp.CN.input$orig.ident==sample.names[i],]
      nearest <- nn2(data=sub.tmp.CN.input[,c("x","y")], k=cell.window)
      for (j in 1 : dim(nearest$nn.idx)[1]){
        if (check.dist) {
          freq <- table(sub.tmp.CN.input[nearest$nn.idx[j,nearest$nn.dists[j,]*pixel.to.micron <= min.dist],"merge.annot"])
        } else {
          freq <- table(sub.tmp.CN.input[nearest$nn.idx[j,],"merge.annot"])
         }
        sub.tmp.CN.input[j, names(freq)] <- freq
      }
      CN.input <- rbind(CN.input,sub.tmp.CN.input)
    }
    CN.input.kmean<-CN.input[,unique(metadata$merge.annot)]
  # identify the cell.neigh CNs via K-means clustering
    MbatchKm<-MiniBatchKmeans(CN.input.kmean,  
                               clusters = cell.neigh,
                               batch_size = 1024,
                               num_init = 3,
                               max_iters = 100,
                               init_fraction = 1,
                               initializer = "kmeans++",
                               early_stop_iter = 10,
                               verbose = FALSE,
                               CENTROIDS = NULL,
                               tol = 1e-04,
                               tol_optimal_init = 0.3,
                               seed = 1)
  # cell type respective frequencies (enrichment score) within each cell neighborhood (CN)
    clus.centroids <- MbatchKm$centroids
    tissue.fract <- table(metadata$merge.annot)/dim(obj.merged)[2]
    tissue.fract <- tissue.fract[unique(metadata$merge.annot)]
    centroids.plus.fraction <- t(apply(clus.centroids, 1 , function(x) x+tissue.fract))
    scale.centroids <- centroids.plus.fraction/rowSums(centroids.plus.fraction)[1]
    data.to.heatmap <- t(apply(scale.centroids, 1 , function(x) x/tissue.fract))
    rownames(data.to.heatmap) <- c(1:cell.neigh)
    annot.row <- data.frame(c(1:cell.neigh))
    rownames(annot.row) <- c(1:cell.neigh)
    colnames(annot.row) <- c("CN")
    CN <- hue_pal()(cell.neigh)
    names(CN) <- levels(as.factor(rownames(annot.row)))
    ann.colors <- list(CN=CN)
    if (check.dist) {suffix <- "_filter"} else {suffix <- "_NOfilter"}
    pdf(file.path(result.dir,paste0("CN_k",cell.neigh,"_w",cell.window,"_identification",suffix,".pdf")),width=10,height=6, useDingbats = FALSE)
    ## heatmap using a custom scale range ##
    colors<-c(-7,seq(-6,6,by=0.01),7)
    my_palette <- c("blue",colorRampPalette(colors = c("blue", "white", "red"))(n = length(colors)-3), "red")
    pheatmap(log2(data.to.heatmap), scale="none",color = my_palette, breaks = colors, main = "Cell type enrichment in cell neighborhoods (CNs)",
             clustering_distance_cols="euclidean", clustering_method="average", treeheight_col = 0,
             cluster_rows=T, cluster_cols=T, border_color="grey90", labels_row=c(1:cell.neigh),
             annotation_row=annot.row, annotation_colors=ann.colors, annotation_names_row=F, annotation_legend=F)
    print(recordPlot())
  # predict cell neighborhoods (CNs) for each cell
    cell.CN <- predict_MBatchKMeans(CN.input.kmean, MbatchKm$centroids, fuzzy = FALSE)
    obj.merged <- AddMetaData(obj.merged, cell.CN, col.name = paste0("neighborhood",cell.neigh))
  # plot for each patient the percent of total cells allocated to each neighborhood
    CN.counts <- table(obj.merged[[]]$orig.ident,obj.merged[[]][,paste0("neighborhood",cell.neigh)])
    CN.freq <- CN.counts/rowSums(CN.counts)
    data.to.plot <- data.frame(CN.freq[sample.names,])
    colnames(data.to.plot) <- c("Sample","CN","Frequency")
    p2 <- ggplot(data.to.plot,aes(x=CN,y=Frequency,color=CN,fill=CN))+
      geom_point(position=position_jitterdodge(jitter.width=1,dodge.width=1),size=2.5, alpha=0.45)+
      stat_summary(fun.data = "mean_cl_boot", colour = "blue", size = 0.6, alpha = 0.75)+
      scale_y_continuous(breaks=seq(0,1,0.1))+
      labs(x="Cell Neighborhood",y="Frequency of Cell Neighborhood")+
      ggtitle("Frequencies of CNs in patients")+
      theme_light()+
      theme(panel.grid.minor = element_blank())+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.4))+
      theme(legend.position="none")+
      geom_text_repel(data=data.to.plot,
                      aes(label=Sample),
                      colour = 'grey50',
                      size=2.5,
                      segment.color="grey50", 
                      segment.size=0.07,
                      min.segment.length = 0.2,
                      max.overlaps=15)
    print(p2)
    # pie charts of CN content in cell types and distribution of cell types in CNs
    palette = c("#fde725", "#21918c", "#440154")
    size.range = c(0,10)
    rotate.x.text = FALSE
    composition<-table(cbind(obj.merged[["merge.annot"]],obj.merged[[paste0("neighborhood",cell.neigh)]]))
    cells.in.CNs <- as.data.frame(composition/rowSums(composition))
    CN.composition <- as.data.frame(t(t(composition)/colSums(composition)))
    plot.out <- data.frame("rows" = CN.composition[,paste0("neighborhood",cell.neigh)],
                     "cols" = CN.composition[,"merge.annot"],
                     "cells.in.CNs" = cells.in.CNs$Freq,
                     "CN.composition" = CN.composition$Freq) 
    p3 <- ggplot(plot.out,aes(x=cols, y=rows, size=CN.composition, fill=cells.in.CNs)) +
      geom_point(shape=21, color = "black") + 
      scale_fill_gradient2(high = palette[1], mid = palette[2], low = palette[3],midpoint = 0.5) + 
      scale_size_continuous(range = size.range) + 
      xlab("Cell type") + 
      ylab("Cellular neighborhood") +
      # theme_minimal_grid() +
      theme_light()+
      theme(axis.text.x = element_text(size = 10, angle = ifelse(rotate.x.text == TRUE, 90, 0))) + 
      labs(size = "Fraction in CN", fill = "Type content") + 
      coord_flip() + 
      ggtitle(paste("Cellular neighborhoods", " vs ", "Cell types"))
    print(p3)
    dev.off()

###########################################################
### Save Seurat object with CN analysis and compositions
###########################################################
  saveRDS(obj.merged, file = file.path(result.dir,paste0("CN_k",cell.neigh,"_w",cell.window,suffix,".rds"))) 
  out.CN <- data.frame(metadata[,keep.cols],CN.input.kmean,cell.CN)
  colnames(out.CN)<- c(keep.cols,unique(metadata$merge.annot),paste0("neighborhood",cell.neigh))
  write.table(out.CN,file.path(result.dir,paste0("CN_k",cell.neigh,"_w",cell.window,suffix,".csv")),sep=',',row.names = F)
  header_st <- createStyle(textDecoration = "Bold")
  wb<-createWorkbook(title="results")
  modifyBaseFont(wb, fontSize = 14, fontColour = "black", fontName = "Arial")
  addWorksheet(wb,paste0(cell.window, "neighbors; ",cell.neigh," CNs"))
  writeData(wb, sheet = paste0(cell.window, "neighbors; ",cell.neigh," CNs"), out.CN, rowNames = F,headerStyle = header_st)
  addWorksheet(wb,paste0("frequencies in patients"))
  writeData(wb, sheet = paste0("frequencies in patients"), CN.freq[sample.names,], rowNames = T,headerStyle = header_st)
  saveWorkbook(wb, file.path(result.dir,paste0("CN_k",cell.neigh,"_w",cell.window,suffix,".xlsx")), overwrite = TRUE)  
    
###########################################################
### draw the cell neighborhoods on the slides for samples 
###########################################################
  if (cell.neigh <=12){
    cols<-brewer.pal(cell.neigh, "Paired")
  } else {
    getPalette <-colorRampPalette(brewer.pal(12, "Paired"))
    cols <- getPalette(cell.neigh)
  }
  names(cols)<-c(1:cell.neigh)
  pdf(file.path(result.dir,paste0("CN_k",cell.neigh,"_w",cell.window,"_slides",suffix,".pdf")),width=12, height=10, useDingbats=T)
  for (i in 1: length(sample.names)){
    sample.dir <- file.path(main.path,"results","single samples",sample.names[i])
    tmp.obj<-readRDS(file.path(sample.dir,paste0(sample.names[i],"_annoted.rds")))
    tmp.obj@meta.data$orig.ident<-sample.names[i]
    tmp.obj <- subset(x = tmp.obj, subset = merge.annot !="artifact")
    sample.CN <- as.factor(out.CN[out.CN$orig.ident==sample.names[i],paste0("neighborhood",cell.neigh)])
    tmp.obj <- AddMetaData(tmp.obj, sample.CN, col.name = paste0("neighborhood",cell.neigh))
    tmp.col <- cols[levels(tmp.obj@meta.data[,paste0("neighborhood",cell.neigh)])]
    print(ImageDimPlot(tmp.obj, fov = sample.names[i], group.by = paste0("neighborhood",cell.neigh),
                       border.color = "NA",cols = tmp.col,size = 1) +
            NoGrid() + 
            scale_y_reverse() +
            ggtitle(sample.names[i]))
    print(ImageDimPlot(tmp.obj, fov = sample.names[i], group.by = paste0("neighborhood",cell.neigh),
                       split.by = paste0("neighborhood",cell.neigh), border.color = "NA",cols = tmp.col,size = 1) +
            NoGrid() + 
            scale_y_reverse() +
            ggtitle(sample.names[i]))
  }
  dev.off()

     