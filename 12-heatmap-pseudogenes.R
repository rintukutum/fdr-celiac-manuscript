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

