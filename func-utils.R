getDetailGeneInfo <- function(entrezIDs){
  entrezIDs <- unique(entrezIDs)
  q.entrez <- na.omit(entrezIDs)
  s.entrez <- q.entrez[which(!duplicated(q.entrez))]
  library(rentrez)
  parts <- seq(1,length(s.entrez),by = 200)
  e_summary <- list()
  cat(paste('No. of chunks\t',length(parts),'\n',sep=''))
  for(i in 1:length(parts)){
    if(i != length(parts)){
      e_summary[[i]] <- entrez_summary(db='gene',
                                       id=as.character(s.entrez[parts[i]:(parts[i+1]-1)]))
    }else{
      e_summary[[i]] <- entrez_summary(db='gene',
                                       id=as.character(s.entrez[parts[i]:length(s.entrez)]))
    }
    cat(paste('chunk\t',i,'\n',sep=''))
  }
  #------
  e_su <- list()
  m <- 1
  for(i in 1:length(e_summary)){
    if(i == 1){
      n <- 200  
    }else{
      n <- n+length(e_summary[[i]])
    }
    e_su[m:n] <- e_summary[[i]]
    m <- n+1
  }
  #-----
  e_header <- c('uid','name','description','status',
                'currentid','chromosome','otheraliases',
                'nomenclaturesymbol','otherdesignations',
                'nomenclaturename','summary','organism'
  )
  e_dat <- lapply(e_su,function(x){
    unlist(x[e_header])
  })
  outHeader <- unique(unlist(lapply(e_dat,names)))
  outDF <- data.frame(matrix(data = NA,
                             ncol = length(outHeader),
                             nrow = length(e_dat)),
                      stringsAsFactors = FALSE)
  colnames(outDF) <- outHeader
  for(i in 1:length(e_dat)){
    outDF[i,names(e_dat[[i]])] <- e_dat[[i]]
  }
  rownames(outDF) <- outDF$uid
  #-----
  # df.entrez <- data.frame(probeID = names(entrezIDs),
  #                         entrezID = entrezIDs,
  #                         stringsAsFactors = FALSE)
  # 
  # library(plyr)
  # t__ <-ddply(.data = df.entrez,
  #             .variables = 'probeID',
  #             .fun = function(x){
  #               df_o <- outDF[as.character(x$entrezID),]
  #               df_o <- cbind(x,df_o)
  #             })
  
  return(outDF)
}
#-----
# Boruta
smartBoruta <- function(x,y,pValue=0.01,maxRuns=100){
  library(Boruta)
  # x = train.dat[,-1]
  # y = train.dat[,1]
  #----
  # split x into bin size f 100
  bin100.start <- seq(1,ncol(x),by = 100)
  bin100.end <- c(bin100.start[-1]-1,
                  ncol(x))
  run.bins <-list()
  for(i in 1:length(bin100.start)){
    run.bins[[i]] <- c(start=bin100.start[i],
                       end=bin100.end[[i]])
  }
  #----
  bt.finalDecision <- lapply(run.bins, function(bin){
    x.dat <- x[,c(bin['start']:bin['end'])]
    bt <- Boruta(x=x.dat,y=y,pValue=pValue,maxRuns=maxRuns)
    return(bt$finalDecision)
  })
  return(bt.finalDecision)
}
getRowColOrder <- function(
  mat
){
  order_r <- hclust(dist(mat))$order
  order_c <- hclust(dist(t(mat)))$order
  return(list(
    r = order_r,
    c = order_c
  ))
}
####################
##### generate heatmap
generateHeatMap <- function(
  normal_probe_data,
  annot.df,
  gene.symbol,
  layout_c,
  dendo = 'column'
){
  idx.gene <- which(annot.df$name %in% gene.symbol)
  probes <- rownames(annot.df)[idx.gene]
  genes <- annot.df[probes,'name']
  mat <- normal_probe_data$df[probes,]
  dat.zScore <- t(apply(mat,1,scale))
  colnames(dat.zScore) <- normal_probe_data$targets
  rownames(dat.zScore) <- genes
  mat2 <- dat.zScore
  orderRC <- getRowColOrder(mat2)
  mat3 <- mat2[orderRC$r,orderRC$c]
  library(dendextend)
  cols <- c(
    c(
      'celiac' = '#217844ff', #green
      'control'= '#0088aaff', #blue, 
      'FDR' = '#aa4400ff') # brown, 
  )[unique(colnames(mat3))]
  
  pchs <- c(
    'celiac' = 0, # CeD
    'control'= 1, # control
    'FDR'= 2 # FDR
  )[unique(colnames(mat3))]
  Colv  <- t(mat3) %>% 
    dist %>% 
    hclust %>% 
    as.dendrogram %>% 
    set("branches_lwd", 1) %>%
    set("leaves_pch", rep(pchs, each=12)) %>%
    set("leaves_cex", 1.5) %>%
    set("leaves_col", rep(cols, each=12))
  ColV <- color_branches(
    dend = Colv,
    k=3,
    col = cols
  )
  cols.gentleman <- function() {
    library(RColorBrewer)
    hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
    return(rev(hmcol))
  }
  #layout.show(layout(mat = layout_c$lmat, widths = layout_c$lwid, heights = layout_c$lhei))
  gplots::heatmap.2(
    mat3,
    offsetCol = 0.25,
    col=cols.gentleman(),
    lmat = layout_c$lmat,
    lwid = layout_c$lwid,
    lhei = layout_c$lhei,
    Colv = ColV,
    dendrogram = dendo,
    tracecol = NULL,
    trace = "none",density.info = 'none',
    key.title ='z-score',
    offsetRow = -0.3,
    key.xlab = 'z-score'
  )
}
layoutConfig <- function(
  vals = c(
    c(0,3),
    c(0,5),
    c(2,1),
    c(0,4)),
  nrow = 4,
  ncol = 2,
  lwid = c(
    0.18,
    3.0
  ),
  lhei = c(
    0.4,
    0.025,
    5,
    0.7
  )
){
  mat <- matrix(
    vals,
    nrow = nrow,
    ncol = ncol,
    byrow = TRUE
  )
  layout_c <- list(
    lmat = mat,
    lwid = lwid,
    lhei = lhei
  )
  return(layout_c)
}
##### generate heatmap
generateHeatMapPair <- function(
  normal_probe_data,
  probe_ID,
  idx,
  layout_c
){
  mat <- normal_probe_data$df[probe_ID,idx]
  dat.zScore <- t(apply(mat,1,scale))
  colnames(dat.zScore) <- normal_probe_data$targets[idx]
  mat <- dat.zScore
  orderRC <- getRowColOrder(mat)
  mat <- mat[orderRC$r,orderRC$c]
  library(dendextend)
  cols <- c(
    c(
      'celiac' = '#217844ff', #green
      'control'= '#0088aaff', #blue, 
      'FDR' = '#aa4400ff') # brown, 
  )[unique(colnames(mat))]
  
  pchs <- c(
    'celiac' = 0, # CeD
    'control'= 1, # control
    'FDR'= 2 # FDR
  )[unique(colnames(mat))]
  Colv  <- t(mat) %>% 
    dist %>% 
    hclust %>% 
    as.dendrogram %>% 
    set("branches_lwd", 1) %>%
    set("leaves_pch", rep(pchs, each=12)) %>%
    set("leaves_cex", 1.5) %>%
    set("leaves_col", rep(cols, each=12))
  ColV <- color_branches(
    dend = Colv,
    k=2,
    col = cols
  )
  cols.gentleman <- function() {
    library(RColorBrewer)
    hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
    return(rev(hmcol))
  }
  #layout.show(layout(mat = layout_c$lmat, widths = layout_c$lwid, heights = layout_c$lhei))
  gplots::heatmap.2(
    mat,
    offsetCol = 0.25,
    col=cols.gentleman(),
    lmat = layout_c$lmat,
    lwid = layout_c$lwid,
    lhei = layout_c$lhei,
    Colv = ColV,
    dendrogram = 'column',
    tracecol = NULL,
    trace = "none",density.info = 'none',
    key.title ='z-score',
    labRow = '',
    offsetRow = -0.3,
    key.xlab = 'z-score'
  )
}
########################
################# GO
getEntrezSig <- function(
  x,
  type = NULL
){
  if(is.null(type)){
    stop('provide up or down')
  }else{
    if(tolower(type) == 'up'){
      idx.sig <- x$log2.FC >= 1
    }
    if(tolower(type) == 'down'){
      idx.sig <- x$log2.FC <= -1
    }
    return(x$EntrezID[idx.sig])
  }
}

doEnrichGO <- function(
  x,
  pval.cutoff = 0.05,
  qval.cutoff = 0.2){
  x <- as.character(x)
  library(clusterProfiler)
  go_bp <- enrichGO(
    gene = x,
    OrgDb = 'org.Hs.eg.db',
    ont = 'BP',
    minGSSize = 1,
    qvalueCutoff = qval.cutoff,
    pvalueCutoff = pval.cutoff,
    readable = TRUE
  )
  go_bp <- go_bp@result
  go_cc <- enrichGO(
    gene = x,
    OrgDb = 'org.Hs.eg.db',
    ont = 'CC',
    minGSSize = 1,
    readable = TRUE,
    qvalueCutoff = qval.cutoff,
    pvalueCutoff = pval.cutoff
  )
  go_cc <- go_cc@result
  go_mf <- enrichGO(
    gene = x,
    OrgDb = 'org.Hs.eg.db',
    ont = 'MF',
    minGSSize = 1,
    qvalueCutoff = qval.cutoff,
    pvalueCutoff = pval.cutoff,
    readable = TRUE
  )
  go_mf <- go_mf@result
  return(
    list(
      'BP' = go_bp,
      'CC' = go_cc,
      'MF' = go_mf
    )
  )
}