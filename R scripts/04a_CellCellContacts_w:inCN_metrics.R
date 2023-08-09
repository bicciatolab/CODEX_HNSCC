####################################################################
### Script to compute differential pairwise cell-cell contacts 
### (cell interaction analysis) within 2 different CNs
####################################################################

library(Seurat)
library(RTriangle)
library(data.table)
library(dplyr)
library(REAT)
library(ggplot2)
library(ggrepel)
library(openxlsx)

project <- "HNSCC"

### main settings ###
  sample.name <- "s2" #"s2","s3","s4","s6","s7","s10","s11","s12","s13"
  # do edit "folder1","folder2","..." with your folder path
  main.path <- file.path("folder1","folder2","...",project)
  self <- "Self" # either "Self" or "noSelf" to include or exclude self-self contacts
  cell.neigh <- 14   # number of total cell neighbors CN c(5,7,10)
  cell.window <- 10  # number of neighboring cells c(5,7,10)
  CN.to.test<-c("TLS1","TLS2") #"Hot Tumor","Cold Tumor"
  
# sub-directory creation #
  data.dir <- file.path(main.path,"results","neighborhood analysis",paste0("CN",cell.neigh,"_window",cell.window))
  result.dir <- file.path(main.path,"results","cell-cell contacts",
                          paste0("within CN (",paste0("k",cell.neigh,"_w",cell.window),")"),
                          "diff interactions")
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
### Identification of cell neighbors of the first (immediate) tier of proximity
### and quantification of llh, relfreq and CLQ in both CNs for every sample   
######################################################################################
  # construct the file for cell contact quantification
  metadata <- data.frame(rownames(obj.merged[[]]),obj.merged[[]])
  colnames(metadata)[1] <- "EventID"
  metadata <- metadata[,c("EventID","orig.ident","x","y","merge.annot",paste0("CN",cell.neigh,"_merged"))]
  # check samples in which there exist both CNs to be compared
  samples.to.use <- intersect(unique(metadata[metadata[,paste0("CN",cell.neigh,"_merged")]==CN.to.test[1],"orig.ident"]),
                              unique(metadata[metadata[,paste0("CN",cell.neigh,"_merged")]==CN.to.test[2],"orig.ident"]))
  llh.freq.CLQ <- NULL
  for (i in 1:length(samples.to.use)){
    sample.llh.freq.CLQ <- NULL
    for (k in 1:length(CN.to.test)){
      sample.data <- metadata[metadata$orig.ident==samples.to.use[i]&metadata[,paste0("CN",cell.neigh,"_merged")]==CN.to.test[k],]
      # check if two cells share the same x/y coordinates
      num.duplicates <- dim(sample.data[duplicated(sample.data[,c("x","y")]),])[1]
      if (num.duplicates==0){print(paste0("No cells sharing the same x/y coordinates in ",CN.to.test[k], " of ",samples.to.use[i]))} else
      {print(paste0(as.character(num.duplicates)," cells sharing the same x/y coordinates in ", CN.to.test[k], " of ", samples.to.use[i]))}
      sample.RT  <- triangulate(pslg(sample.data[,c("x","y")]), Y = T, D = T)
      # Delauney graph calculation using RTriangle package. Duplicate A-B contacts to get also B-A contacts
      sample.resultRT <- rbind(sample.RT$E,sample.RT$E[,c(2,1)])
      sample.resultRT <- data.frame(sample.data[sample.resultRT[,1],"x"],
                                    sample.data[sample.resultRT[,1],"y"],
                                    sample.data[sample.resultRT[,2],"x"],
                                    sample.data[sample.resultRT[,2],"y"],
                                    sample.data[sample.resultRT[,1],"merge.annot"],
                                    sample.data[sample.resultRT[,2],"merge.annot"],
                                    rep(sample.names[i],dim(sample.resultRT)[1]),
                                    rep(CN.to.test[k],dim(sample.resultRT)[1]))
      colnames(sample.resultRT) <- c("x1","y1","x2","y2","cell1type","cell2type","samplename","CN")
      ## calculate contact likelihood ratio and relative frequencies
      sample.resultRT <- as.data.table(sample.resultRT)
      # N_12: interaction counts for cell_Type_1 with cell_Type_2
      cn.llh.freq.CLQ <- sample.resultRT[, .(N_12 = .N), .(cell1type, cell2type, CN)] %>%
        # N_1: interaction counts for cell_Type_1
        left_join(sample.resultRT[, .(N_1 = .N), .(cell1type, CN)],by = c('cell1type', 'CN')) %>% 
        # N_2: interaction counts for cell_Type_2
        left_join(sample.resultRT[, .(N_2 = .N), .(cell2type, CN)],by = c('cell2type', 'CN')) %>% 
        # N_total: total interaction counts
        left_join(sample.resultRT[, .(N_total = .N), CN], by = 'CN') %>%
        # likelihood ratio for interactions (N_12*N_total)/(N_1*N_2)
        mutate(llh = (1.0 * N_12 * N_total) / (1.0 * N_1 * N_2)) %>%
        # relative interaction frequency N_12/N_1
        mutate(relfreq = 1.0 * N_12/N_1) %>%
        as.data.table
      ## add sample name and interactions
      cn.llh.freq.CLQ <- cn.llh.freq.CLQ %>%
        mutate(sample = rep(samples.to.use[i],dim(cn.llh.freq.CLQ)[1])) %>%
        mutate(inter = paste0(cn.llh.freq.CLQ$cell1type,"-->",cn.llh.freq.CLQ$cell2type)) %>%
        as.data.table
      sample.llh.freq.CLQ <- rbind(sample.llh.freq.CLQ,cn.llh.freq.CLQ)
    }
    llh.freq.CLQ <- rbind(llh.freq.CLQ,sample.llh.freq.CLQ)
  }

########################################################################################
### Quantification of interaction LLH, Relative Frequency and Co-localization Quotient 
### of cell-cell interactions present in at least n samples per each CNs
########################################################################################
  n.samples <- 2
  # cell-cell interactions in each of the 2 CNs
  CN.1 <- llh.freq.CLQ[llh.freq.CLQ$CN==CN.to.test[1],]
  CN.2 <- llh.freq.CLQ[llh.freq.CLQ$CN==CN.to.test[2],]
  # cell-cell interactions that are present in at least n samples in each of the two CNs
  inter.1 <- names(table(CN.1$inter)[table(CN.1$inter) >= n.samples])
  inter.2 <- names(table(CN.2$inter)[table(CN.2$inter) >= n.samples])
  # common cell-cell interactions in at least n samples in each of the two CNs
  common.inter <- intersect(inter.1,inter.2)
  CN.1.inter <- setdiff(CN.1$inter,CN.2$inter) # interactions present only in the first CN
  CN.2.inter <- setdiff(CN.2$inter,CN.1$inter) # interactions present only in the second CN
  diff.llh <- matrix(data = NA, nrow = length(common.inter), ncol = 4)
  diff.freq <- diff.llh
  diff.clq <- diff.llh
  for (j in 1:length(common.inter)){
    diff.llh[j,1] <- mean(CN.1[inter %in% common.inter[j],llh])/mean(CN.2[inter %in% common.inter[j],llh])
    diff.llh[j,2] <- log2(diff.llh[j,1])
    diff.llh[j,3] <- t.test(CN.1[inter %in% common.inter[j],llh],CN.2[inter %in% common.inter[j],llh])$p.value
    diff.freq[j,1] <- mean(CN.1[inter %in% common.inter[j],relfreq])/mean(CN.2[inter %in% common.inter[j],relfreq])
    diff.freq[j,2] <- log2(diff.freq[j,1])
    diff.freq[j,3] <- t.test(CN.1[inter %in% common.inter[j],relfreq],CN.2[inter %in% common.inter[j],relfreq])$p.value
    diff.clq[j,1] <- mean(CN.1[inter %in% common.inter[j],CLQ])/mean(CN.2[inter %in% common.inter[j],CLQ])
    diff.clq[j,2] <- log2(diff.clq[j,1])
    diff.clq[j,3] <- t.test(CN.1[inter %in% common.inter[j],CLQ],CN.2[inter %in% common.inter[j],CLQ])$p.value
  }
  rownames(diff.llh) <- common.inter
  rownames(diff.freq) <- common.inter
  rownames(diff.clq) <- common.inter
  colnames(diff.llh) <- c("FC","log2FC","p.value","p.adj")
  colnames(diff.freq) <- c("FC","log2FC","p.value","p.adj")
  colnames(diff.clq) <- c("FC","log2FC","p.value","p.adj")
  diff.llh <- as.data.frame(diff.llh)
  diff.freq <- as.data.frame(diff.freq)
  diff.clq <- as.data.frame(diff.clq)
  diff.llh$p.adj <- p.adjust(diff.llh$p.value, method = "BH")
  diff.freq$p.adj <- p.adjust(diff.freq$p.value, method = "BH")
  diff.clq$p.adj <- p.adjust(diff.clq$p.value, method = "BH")
  diff.llh$FC<-ifelse(diff.llh$FC>1,diff.llh$FC,-1/diff.llh$FC)
  diff.freq$FC<-ifelse(diff.freq$FC>1,diff.freq$FC,-1/diff.freq$FC)
  diff.clq$FC<-ifelse(diff.clq$FC>1,diff.clq$FC,-1/diff.clq$FC)
  
########################################################################################
### Volcano plot of LLH, Relative Frequency and Co-localization Quotient 
### of cell-cell interactions present in at least n samples per each CNs
########################################################################################
  metric.to.plot <- "RelFreq"
  if (metric.to.plot == "LLH"){volc.table <- diff.llh}
  if (metric.to.plot == "RelFreq"){volc.table <- diff.freq}
  if (metric.to.plot == "CLQ"){volc.table <- diff.clq}
  fold.change <- 2
  FDR.cutoff <- 0.05
  FC.cutoff <- log2(fold.change)
  # contruct the matrix for the volcano plot
  volc.table$Significance <- "NS"
  volc.table$Significance[volc.table$p.adj<=FDR.cutoff] <- "FDRonly"
  volc.table$Significance[volc.table$p.adj<=FDR.cutoff&volc.table$log2FC>=FC.cutoff] <- "up"
  volc.table$Significance[volc.table$p.adj<=FDR.cutoff&volc.table$log2FC<=-FC.cutoff] <- "down"
  table(volc.table$Significance)
  volc.table$Significance <- factor(volc.table$Significance, levels=c("NS", "FDRonly", "up","down"))
  # volcano plot 
  pdf(file.path(result.dir,paste0(CN.to.test[1],".vs.",CN.to.test[2],"_",metric.to.plot,".pdf")), height=8, width=12, useDingbats=FALSE)
  p1 <- ggplot(volc.table, aes(x=log2FC, y=-log10(p.adj))) +
    geom_point(aes(color=Significance)) +
    scale_color_manual(values=c(NS="gray80",FDRonly="gray80",up="red",down="blue"),breaks=c("up","down"))+
    #xlim(-3, 3) +
    #ylim(0, 9) +
    xlab(bquote(~log["2"]~"(FoldChange)"))+
    ylab(bquote(~-log["10"]~"(P value)")) +
    ggtitle(paste0(CN.to.test[1]," vs. ",CN.to.test[2]," (",metric.to.plot,")" )) +
    theme_light()+
    theme(panel.grid.minor = element_blank())+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.4))+
    theme(legend.position="none")+
    theme(plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
          axis.title=element_text(size=10),
          aspect.ratio=0.5) +
    geom_text_repel(data=volc.table[volc.table$Significance=="up"|volc.table$Significance=="down",],
                    aes(label=rownames(volc.table[volc.table$Significance=="up"|volc.table$Significance=="down",])),
                    size=2.25,
                    segment.color="black", 
                    segment.size=0.07,
                    #direction = "y",
                    hjust = 1,
                    vjust= 0.4,
                    min.segment.length = 0.1,
                    max.overlaps=20)
  print(p1)
  dev.off()
  
########################################################################################
### Write metrics and metrics differentials
########################################################################################
  header_st <- createStyle(textDecoration = "Bold")
  cn.metrics <-createWorkbook(title="interaction metrics in CNs")
  modifyBaseFont(cn.metrics, fontSize = 14, fontColour = "black", fontName = "Arial")
  addWorksheet(cn.metrics,"metrics")
  llh.freq.CLQ[,c("llh","relfreq","CLQ")]<-round(llh.freq.CLQ[,c("llh","relfreq","CLQ")],3)
  writeData(cn.metrics, sheet = "metrics", llh.freq.CLQ, rowNames = F,headerStyle = header_st)
  addWorksheet(cn.metrics,"diff on LLH")
  diff.llh.out<-data.frame("interaction" = rownames(diff.llh),
                       diff.llh)
  diff.llh.out[,c("FC","log2FC")]<-round(diff.llh.out[,c("FC","log2FC")],2)
  diff.llh.out[,c("p.value","p.adj")]<-round(diff.llh.out[,c("p.value","p.adj")],4)
  writeData(cn.metrics, sheet = "diff on LLH", diff.llh.out, rowNames = F,headerStyle = header_st)
  addWorksheet(cn.metrics,"diff on RelFreq")
  diff.freq.out<-data.frame("interaction" = rownames(diff.freq),
                        diff.freq)
  diff.freq.out[,c("FC","log2FC")]<-round(diff.freq.out[,c("FC","log2FC")],2)
  diff.freq.out[,c("p.value","p.adj")]<-round(diff.freq.out[,c("p.value","p.adj")],4)
  writeData(cn.metrics, sheet = "diff on RelFreq", diff.freq.out, rowNames = F,headerStyle = header_st)
  addWorksheet(cn.metrics,"diff on CLQ")
  diff.clq.out<-data.frame("interaction" = rownames(diff.clq),
                       diff.clq)
  diff.clq.out[,c("FC","log2FC")]<-round(diff.clq.out[,c("FC","log2FC")],2)
  diff.clq.out[,c("p.value","p.adj")]<-round(diff.clq.out[,c("p.value","p.adj")],4)
  writeData(cn.metrics, sheet = "diff on CLQ", diff.clq.out, rowNames = F,headerStyle = header_st)
  saveWorkbook(cn.metrics, file.path(result.dir,paste0(CN.to.test[1],".vs.",CN.to.test[2],"_metrics.xlsx")), overwrite = TRUE)  
  
  
