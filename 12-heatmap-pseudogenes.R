rm(list=ls())
de_probes <- read.csv(
  './data/DE_probes_annotated_pvals_fold_change.csv',
  stringsAsFactors = FALSE,
  row.names = 1)
idx.pseudogenes <- grep('pseudogene', de_probes$description)
de_probes.pseudo <- de_probes[idx.pseudogenes,]
de_probes.gene.desc <- de_probes.pseudo[,c('name','description')]
de_probes.gene.desc <- cbind(
  probeID = rownames(de_probes.gene.desc),
  de_probes.gene.desc)
write.csv(de_probes.gene.desc,
          './data/de_probes_pseudo_genes_description.csv',
          row.names = FALSE)
#-------------
de_probes.pseudo <- de_probes.pseudo[order(de_probes.pseudo$name,decreasing = TRUE),]
load('./data/normal_probe_data.RData')
pseudo.mat <- normal_probe_data$df[rownames(de_probes.pseudo),]
###
source('./func-utils.R')
layout_c <- layoutConfig(
  vals = c(c(0,3), c(0,5), c(2,1), c(0,4)),
  nrow = 4,
  ncol = 2,
  lwid = c(1.5, 3.0),
  lhei = c(0.35, 0.03, 6.5, 0.4)
)
pdf('./figures/12-pseudo-genes-heatmap.pdf',
    height = 22,
    width = 7)
generateHeatMap(
  normal_probe_data = normal_probe_data,
  annot.df = de_probes.pseudo,
  gene.symbol = de_probes.pseudo$name,
  layout_c = layout_c,
  dendo = 'both')
dev.off()

#--------------
# annotate pseudogenes and categories into different 
# categories 

# download.file(
#   'http://pseudogene.org/psicube/data/gencode.v10.pseudogene.txt',
#   './data/gencode.v10.pseudogene.txt')
#####
# pseudogenes.annotate <- readLines('./data/gencode.v10.pseudogene.txt')
# getNameVal <- function(x){
#   xx <- strsplit(x,split = ' ')[[1]]
#   y <- xx[2]
#   names(y) <- xx[1]
#   return(y)
# }
# getInfo <- function(x){
#   x <- strsplit(gsub('\"','',strsplit(x,split='\t')[[1]][9]),split='\\;')[[1]]
#   x <- gsub('^[[:space:]]','',x)
#   out.x <- list()
#   for(i in 1:length(x)){
#     out.x[[i]] <- getNameVal(x[i])
#   }
#   out <- unlist(out.x)
#   return(out)
# }
# pseudo.genes <- list()
# for(i in 1:length(pseudogenes.annotate)){
#   pseudo.genes[[i]] <- getInfo(pseudogenes.annotate[i])
# }
# pseudo.genenames <- c()
# for(i in 1:length(pseudo.genes)){
#   p.gene <- pseudo.genes[[i]]['gene_name']
#   if(length(grep('-',p.gene)) != 0){
#     pseudo.genenames[i] <- strsplit(p.gene,split = '\\-')[[1]][1]
#   }else{
#     if(length(grep('\\.',p.gene)) != 0){
#       pseudo.genenames[i] <- strsplit(p.gene,split = '\\.')[[1]][1]
#     }else{
#       pseudo.genenames[i] <- p.gene
#     }
#   }
# }
# idx.pseu.anno <- which(pseudo.genenames %in% de_probes.gene.desc$name)
# pgo <- sapply(pseudo.genes[idx.pseu.anno],function(x){x['ont']})
# pge <- pseudo.genenames[idx.pseu.anno]
# pgoe <- data.frame(ont = pgo,name = pge,stringsAsFactors = FALSE)
# pgoe.ont <- plyr::dlply(pgoe,'ont')
# ################
# source('./func-utils.R')
# curr.layout <- layoutConfig(
#   vals = c(c(0,3), c(0,5), c(2,1), c(0,4)),
#   nrow = 4,
#   ncol = 2,
#   lwid = c(0.18, 3.0),
#   lhei = c(0.8, 0.03, 5, 1)
#   )
# pseudogene.type <- c(
#   'unitary' = 'PGO:0000004',
#   'duplicated' = 'PGO:0000005'
# )
# pgo.name <- names(pgoe.ont)[1]
# pgo.type <- names(pseudogene.type[pseudogene.type %in% pgo.name])
# file.name <- paste('./figures/12a-pseudogene-',pgo.type,'.pdf')
# pdf(file.name, width = 5, height = 10)
# file.name <- paste('./figures/12a-pseudogene-',pgo.type,'.svg')
# svg(file.name, width = 5, height = 10)
# par(mar=c(5.1,4.1,4.1,4.1))
# generateHeatMap(
#   normal_probe_data = normal_probe_data,
#   annot.df = de_probes.gene.desc,
#   gene.symbol = pgoe.ont[[1]]$name,
#   layout_c = curr.layout
#   )
# dev.off()
# file.name <- paste('./figures/12a-pseudogene-',pgo.type,'.png')
# 
# #############
# pgo.name <- names(pgoe.ont)[2]
# pgo.type <- names(pseudogene.type[pseudogene.type %in% pgo.name])
# pdf(
#   paste(
#     './figures/12b-pseudogene-',
#     pgo.type,
#     '.pdf'
#   ),
#   width = 5,
#   height = 7.5
# )
# par(mar=c(5.1,4.1,4.1,4.1))
# curr.layout$lhei[1] <- 1
# curr.layout$lhei[3] <- 4
# generateHeatMap(
#   normal_probe_data = normal_probe_data,
#   annot.df = de_probes.gene.desc,
#   gene.symbol = pgoe.ont[[2]]$name,
#   layout_c = curr.layout
# )
# dev.off()