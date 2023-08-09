########################################################################################
### Script to compute pairwise cell-cell contacts (cell interaction analysis) within CN
########################################################################################

library(Seurat)
library(RTriangle)
library(paletteer)
library(circlize)
library(gplots)
library(pheatmap)
library(openxlsx)

project <- "HNSCC"

### main settings ###
  sample.name <- "s2" #"s2","s3","s4","s6","s7","s10","s11","s12","s13"
  # do edit "folder1","folder2","..." with your folder path
  main.path <- file.path("folder1","folder2","...",project)
  self <- "Self" # either "Self" or "noSelf" to include or exclude self-self contacts
  cell.neigh <- 14   # number of total cell neighbors CN c(5,7,10)
  cell.window <- 10  # number of neighboring cells c(5,7,10)

# sub-directory creation #
  data.dir <- file.path(main.path,"results","neighborhood analysis",paste0("CN",cell.neigh,"_window",cell.window))
  main.result.dir <- file.path(main.path,"results","cell-cell contacts")
  if (!file.exists(main.result.dir)){dir.create(main.result.dir)}
  result.dir <- file.path(main.path,"results","cell-cell contacts",paste0("within CN (",paste0("k",cell.neigh,"_w",cell.window),")"))
  if (!file.exists(result.dir)){dir.create(result.dir)}

################################################################################
### Load Seurat object and CN annotation 
################################################################################
  obj.merged <- readRDS(file.path(data.dir,paste0("CN_k",cell.neigh,"_w",cell.window,"_NOfilter.rds")))
  CN.annot <- read.table(file.path(data.dir,paste0("CN_k",cell.neigh,"_w",cell.window,"_NOfilter_annot.txt")),
                           sep='\t',header = T,stringsAsFactors = F)
  # CN re-labeling #
  CN.annot.col <- paste0("CN",cell.neigh,"_merged")
  obj.merged[[CN.annot.col]]<-obj.merged[[paste0("neighborhood",cell.neigh)]]
  codex.merge.annot<-obj.merged[[CN.annot.col]]
  for (i in 1:dim(CN.annot)[1]) {
    codex.merge.annot[,CN.annot.col][codex.merge.annot[,CN.annot.col]==CN.annot[i,1]] <- CN.annot[i,CN.annot.col]
  }
  obj.merged[[CN.annot.col]]<-codex.merge.annot[,CN.annot.col]
  
######################################################################################
### Identification of cell neighbors of the first (immediate) tier of proximity and
### quantification of odds ratio of co-occurrence of cell type A and cell type B   
######################################################################################
  self <- "Self" # either "Self" or "noSelf" to include or exclude self-self contacts
  # construct the file for cell contact quantification
  metadata <- data.frame(rownames(obj.merged[[]]),obj.merged[[]])
  colnames(metadata)[1] <- "EventID"
  metadata <- metadata[,c("EventID","orig.ident","x","y",CN.annot.col,"merge.annot")]
  # Delauney graph calculations and log likelihood ratio
  for (j in 1:length(sample.names)){
    sample.data.single <- metadata[metadata$orig.ident==sample.names[j],]
    chord.list <- list()
    heat.list <- list()
    header_st <- createStyle(textDecoration = "Bold")
    cc.concts <-createWorkbook(title="cell-cell contacts")
    cc.llhratios <-createWorkbook(title="likelihodd ratios")
    modifyBaseFont(cc.concts, fontSize = 14, fontColour = "black", fontName = "Arial")
    modifyBaseFont(cc.llhratios, fontSize = 14, fontColour = "black", fontName = "Arial")
    CN.annot.sample <- unique(sample.data.single[,CN.annot.col])
    CN.annot.sample <- unique(CN.annot[,CN.annot.col][CN.annot[,CN.annot.col]%in%CN.annot.sample])
    CN.number.sample <- CN.annot[,1][CN.annot[,CN.annot.col]%in%CN.annot.sample]
    for (i in 1:length(CN.annot.sample)){
      print(paste0(sample.names[j]," - ",CN.annot.sample[i]))
      sample.data <- sample.data.single[sample.data.single[,CN.annot.col]==CN.annot.sample[i],]
      # check if two cells share the same x/y coordinates
      num.duplicates <- dim(sample.data[duplicated(sample.data[,c("x","y")]),])[1]
      if (num.duplicates==0){print(paste0("No cells sharing the same x/y coordinates"))} else
      { print(paste0(as.character(num.duplicates)," cells sharing the same x/y coordinates"))
        rows.to.keep <- setdiff(rownames(sample.data),rownames(sample.data[duplicated(sample.data[,c("x","y")]),]))
        sample.data <- sample.data[rows.to.keep,]
      }
      # Delauney graph calculation using RTriangle package
        sample.RT  <- triangulate(pslg(sample.data[,c("x","y")]), Y = T, D = T)
        # duplicate A-B contacts to get also B-A contacts
        sample.resultRT <- rbind(sample.RT$E,sample.RT$E[,c(2,1)])
        sample.resultRT <- data.frame(sample.data[sample.resultRT[,1],"merge.annot"],sample.data[sample.resultRT[,2],"merge.annot"])
        colnames(sample.resultRT) <- c("cell1type","cell2type")
      # contact counts  
        pivot.ori <- table(sample.resultRT)
        pivot.norm <- t(pivot.ori/rowSums(pivot.ori))
      # chord diagram of contact counts
        data.to.chord <- pivot.ori
        miss.marker <- setdiff(unique(metadata$merge.annot),colnames(data.to.chord))
        if (length(miss.marker) >0){
          miss.marker.mat1 <- matrix(data=0.3,nrow=length(miss.marker),ncol=dim(data.to.chord)[2])
          miss.marker.mat2 <- matrix(data=0.3,nrow=length(unique(metadata$merge.annot)),ncol=length(miss.marker))
          rownames(miss.marker.mat1) <- miss.marker
          colnames(miss.marker.mat2) <- miss.marker
          data.to.chord <- rbind(data.to.chord,miss.marker.mat1)
          data.to.chord <- cbind(data.to.chord,miss.marker.mat2)
          data.to.chord <- as.data.frame(data.to.chord)
          data.to.chord <- data.to.chord[order(row.names(data.to.chord)),]
          data.to.chord <- data.to.chord[,order(names(data.to.chord))]
        }
        if (self == "noSelf"){
          data.to.chord[lower.tri(data.to.chord,diag=T)] <- 0
        } else {
          data.to.chord[lower.tri(data.to.chord,diag=F)] <- 0
        }
        data.to.chord <- as.matrix(data.to.chord)
        # set the colors
        getPalette <- colorRampPalette(paletteer_d("basetheme::dark"))
        grid.col <- getPalette(length(unique(colnames(data.to.chord))))
        color.mat <- t(matrix(data=grid.col,
                              nrow = length(unique(colnames(data.to.chord))),
                              ncol = length(unique(colnames(data.to.chord)))))
        color.mat[] <- Vectorize(adjustcolor)(color.mat, alpha.f = 0.6)
        # chord diagram
        circos.clear()
        par(cex = 0.5, mar = c(0, 0, 4, 0))
        circos.par(gap.degree=4)
        set_track_gap(mm_h(1))
        chordDiagram(data.to.chord,
                     annotationTrack = "grid",
                     preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(data.to.chord))))),
                     grid.col = grid.col,
                     col = color.mat,
                     scale = F)
        circos.track(track.index = 1, panel.fun = function(x, y) {
          circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                      facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                      cex = fontsize(16))
        }, bg.border = NA)
        title(paste0(sample.names[j]," - CN", CN.number.sample[i],": ",CN.annot.sample[i], 
                     " (cell-cell contacts)"), cex.main = 2)
        chord.list[[i]] <- recordPlot()
      # contact likelihood ratio
        summa <- sum(pivot.ori)
        hornorm <- unname(rowSums(pivot.ori))/summa
        lglikelihood <- as.matrix(pivot.ori)
        for (k in 1:dim(pivot.ori)[1]) {
          for (h in 1: dim(pivot.ori)[2]) {
            lglikelihood[k,h] = lglikelihood[k,h]/(hornorm[k]*hornorm[h]*summa)
          }
        }
      # heatmap of contact likelihood ratios  
        lglikelihood.to.heat <-lglikelihood
        miss.marker <- setdiff(unique(metadata$merge.annot),colnames(lglikelihood.to.heat))
        if (length(miss.marker) >0){
          miss.marker.mat1 <- matrix(data=NA,nrow=length(miss.marker),ncol=dim(lglikelihood.to.heat)[2])
          miss.marker.mat2 <- matrix(data=NA,nrow=length(unique(metadata$merge.annot)),ncol=length(miss.marker))
          rownames(miss.marker.mat1) <- miss.marker
          colnames(miss.marker.mat2) <- miss.marker
          lglikelihood.to.heat <- rbind(lglikelihood.to.heat,miss.marker.mat1)
          lglikelihood.to.heat <- cbind(lglikelihood.to.heat,miss.marker.mat2)
          lglikelihood.to.heat <- as.data.frame(lglikelihood.to.heat)
          lglikelihood.to.heat <- lglikelihood.to.heat[order(row.names(lglikelihood.to.heat)),]
          lglikelihood.to.heat <- lglikelihood.to.heat[,order(names(lglikelihood.to.heat))]
        }
        lglikelihood.to.heat[lower.tri(lglikelihood.to.heat)] <- 1
        breaks <- seq(0, 4, length.out = 1000)
        gradient1 <- colorpanel( sum( breaks[-1]<=1 ), "#0030FF", "white" )
        gradient2 <- colorpanel( sum( breaks[-1]>1 ), "white", "#FF5020" )
        hm.colors <- c(gradient1,gradient2)
        pheatmap(lglikelihood.to.heat, scale="none", clustering_method = "average", cluster_rows = F,
                 cluster_cols = F, color = hm.colors, breaks = breaks, border_color = "white",
                 na_col = "grey97", cellwidth=17,cellheight=17,
                 main = paste0(sample.names[j]," - CN", CN.number.sample[i],": ",
                               CN.annot.sample[i], " (likelihood ratios)"))
        heat.list[[i]] <- recordPlot()
      # writing contacts to Excel
        addWorksheet(cc.concts,CN.annot.sample[i])
        cc.excel <- pivot.ori
        miss.marker <- setdiff(unique(metadata$merge.annot),colnames(cc.excel))
        if (length(miss.marker) >0){
          miss.marker.mat1 <- matrix(data='',nrow=length(miss.marker),ncol=dim(cc.excel)[2])
          miss.marker.mat2 <- matrix(data='',nrow=length(unique(metadata$merge.annot)),ncol=length(miss.marker))
          rownames(miss.marker.mat1) <- miss.marker
          colnames(miss.marker.mat2) <- miss.marker
          cc.excel <- rbind(cc.excel,miss.marker.mat1)
          cc.excel <- cbind(cc.excel,miss.marker.mat2)
          cc.excel <- as.data.frame(cc.excel)
          cc.excel <- cc.excel[order(row.names(cc.excel)),]
          cc.excel <- cc.excel[,order(names(cc.excel))]
        }  
        cc.excel[lower.tri(cc.excel,diag=F)] <- ''
        writeData(cc.concts, sheet = CN.annot.sample[i], cc.excel, rowNames = T,headerStyle = header_st)
      # writing likelihood ratios to Excel
        addWorksheet(cc.llhratios,CN.annot.sample[i])
        ll.excel <- round(lglikelihood,2)
        miss.marker <- setdiff(unique(metadata$merge.annot),colnames(ll.excel))
        if (length(miss.marker) >0){
          miss.marker.mat1 <- matrix(data='',nrow=length(miss.marker),ncol=dim(ll.excel)[2])
          miss.marker.mat2 <- matrix(data='',nrow=length(unique(metadata$merge.annot)),ncol=length(miss.marker))
          rownames(miss.marker.mat1) <- miss.marker
          colnames(miss.marker.mat2) <- miss.marker
          ll.excel <- rbind(ll.excel,miss.marker.mat1)
          ll.excel <- cbind(ll.excel,miss.marker.mat2)
          ll.excel <- as.data.frame(ll.excel)
          ll.excel <- ll.excel[order(row.names(ll.excel)),]
          ll.excel <- ll.excel[,order(names(ll.excel))]
        }
        ll.excel[lower.tri(ll.excel,diag=F)] <- ''
        writeData(cc.llhratios, sheet = CN.annot.sample[i], ll.excel, rowNames = T,headerStyle = header_st)
    }
  # saving chord diagrams  
    pdf(file.path(result.dir,paste0(sample.names[j],"_","CellCellContacts_witinCN_",self,".pdf")),
        width=7,height=7, useDingbats = FALSE)
    for (rr in 1:length(CN.annot.sample)){
      print(chord.list[[rr]])
    }
    dev.off()
  # saving heatmaps of likelihood ratios   
    pdf(file.path(result.dir,paste0(sample.names[j],"_","LikelihoodRatios_withinCN.pdf")),
        width=7,height=7, useDingbats = FALSE)
    for (rr in 1:length(CN.annot.sample)){
      print(heat.list[[rr]])
    }
    dev.off()
  # saving cell-cell contacts and likelihood ratios
    saveWorkbook(cc.concts, file.path(result.dir,paste0(sample.names[j],"_","CellCellContacts_withinCN.xlsx")),
                 overwrite = TRUE)  
    saveWorkbook(cc.llhratios, file.path(result.dir,paste0(sample.names[j],"_","LikelihoodRatios_withinCN.xlsx")),
                 overwrite = TRUE)
  } 
  
