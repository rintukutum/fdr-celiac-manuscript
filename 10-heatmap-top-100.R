rm(list=ls())
var_imp <- read.csv(
  './data/rf_variable_importance.csv',
  stringsAsFactors = FALSE,
  row.names = 1)
top100probes <- rownames(var_imp)[order(var_imp$MeanDecreaseAccuracy,
                                        decreasing = TRUE)[1:100]]
anno.de.probes <- read.csv(
  './data/differentially_expressed_probes_annotation_entrez.csv',
  stringsAsFactors = FALSE
)
rownames(anno.de.probes) <- anno.de.probes$probeID
probe2gene <- anno.de.probes[top100probes,'name']
load('./data/normal_probe_data.RData')
source('./func-utils.R')
layout_c <- layoutConfig(
  vals = c(c(0,3), c(0,5), c(2,1), c(0,4)),
  nrow = 4,
  ncol = 2,
  lwid = c(0.18, 3.0),
  lhei = c(0.35, 0.03, 3, 0.4)
)
pdf('./figures/11-top100-probes-random-forest.pdf',
    width = 5,
    height = 12)
generateHeatMap(
  normal_probe_data = normal_probe_data,
  annot.df = anno.de.probes,
  gene.symbol = probe2gene,
  layout_c = layout_c
)
dev.off()
###########################
