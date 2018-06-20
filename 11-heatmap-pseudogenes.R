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

############################## cluster members manual extraction
annot.df <- de_probes.pseudo
gene.symbol <- de_probes.pseudo$name
idx.gene <- which(annot.df$name %in% gene.symbol)
probes <- rownames(annot.df)[idx.gene]
genes <- annot.df[probes,'name']
mat <- normal_probe_data$df[probes,]
rownames(mat) <- paste(probes,'.',genes,sep='')
dat.zScore <- t(apply(mat,1,scale))
ht2 <- heatmap(dat.zScore)
genes_order <- sapply(rev(rownames(mat)[ht2$rowInd]),
	function(x){strsplit(x,split='\\.')[[1]][2]}
)
probes_order <- sapply(rev(rownames(mat)[ht2$rowInd]),
	function(x){strsplit(x,split='\\.')[[1]][1]}
)

end_cluster_1 <- grep('SRRM1P1',genes_order)
end_cluster_2 <- grep('RPL13AP5',genes_order)[2]
end_cluster_3 <- grep('RPS15P4',genes_order)
end_cluster_4 <- length(genes_order)

cluster_order <- rep(NA,length(genes_order))
cluster_order[1:end_cluster_1] <- 'A'
cluster_order[(end_cluster_1+1):end_cluster_2] <- 'B'
cluster_order[(end_cluster_2+1):end_cluster_3] <- 'C'
cluster_order[(end_cluster_3+1):end_cluster_4] <- 'D'

cluster_info <- data.frame(
	probeID = probes_order,
	geneID = genes_order,
	cluster = cluster_order,
	stringsAsFactors = FALSE
)
pseudo_gene_anno <- read.csv(
	'./data/de_probes_pseudo_genes_description.csv',
	stringsAsFactors=FALSE)
rownames(pseudo_gene_anno) <- pseudo_gene_anno$probeID

cluster_info_n <- data.frame(
	cluster_info,
	description = pseudo_gene_anno[cluster_info$probeID,'description'],
	stringsAsFactors=FALSE
)
write.csv(cluster_info_n,
	'./data/de_probes_pseudo_genes_cluster_ID_description.csv',
	row.names=FALSE
)
