rm(list=ls())
load('./data/anno.de.df.RData')
load('./data/normal_probe_data.RData')
de_probes <- read.csv(
  './data/DE_probes_annotated_pvals_fold_change.csv',
  stringsAsFactors = FALSE)
idx.probe <- de_probes$celiac_vs_FDR.pval_adjusted <= 0.05 & abs(de_probes$celiac_vs_FDR.fc) >= 1
idx.FDR <- which(normal_probe_data$targets == "FDR")
idx.CeD <- which(normal_probe_data$targets == "celiac")
idx.Ctrl <- which(normal_probe_data$targets == "control")
##--- FDR vs CeD
idx.probe <- de_probes$celiac_vs_FDR.pval_adjusted <= 0.05 & abs(de_probes$celiac_vs_FDR.fc) >= 1
FDR.CeD.probes <- de_probes$probeID[idx.probe]
# compute z-score
source('./func-utils.R')
layout_c <- layoutConfig(
  vals = c(c(0,3), c(0,5), c(2,1), c(0,4)),
  nrow = 4,
  ncol = 2,
  lwid = c(0.18, 3.0),
  lhei = c(0.35, 0.03, 5, 0.4)
)
pdf('./figures/3a-heatmap_FDR_vs_CeD.pdf',
    width = 6,
    height = 18)
generateHeatMapPair(
  normal_probe_data = normal_probe_data,
  probe_ID = FDR.CeD.probes,
  idx = c(idx.FDR,idx.CeD),
  layout_c = layout_c
  )
dev.off()
##--- Ctrl vs CeD
# compute z-score
idx.probe <- de_probes$celiac_vs_control.pval_adjusted <= 0.05 & abs(de_probes$celiac_vs_control.fc) >= 1
Ctrl.CeD.probes <- de_probes$probeID[idx.probe]
# compute z-score
layout_c <- layoutConfig(
  vals = c(c(0,3), c(0,5), c(2,1), c(0,4)),
  nrow = 4,
  ncol = 2,
  lwid = c(0.18, 3.0),
  lhei = c(0.35, 0.03, 3.5, 0.6)
)
pdf('./figures/3b-heatmap_Ctrl_vs_CeD.pdf',
    width = 6,
    height = 12)
generateHeatMapPair(
  normal_probe_data = normal_probe_data,
  probe_ID = Ctrl.CeD.probes,
  idx = c(idx.Ctrl,idx.CeD),
  layout_c = layout_c
)
dev.off()

##--- Ctrl vs FDR
# compute z-score
idx.probe <- de_probes$FDR_vs_control.pval_adjusted <= 0.05 & abs(de_probes$FDR_vs_control.fc) >= 1
Ctrl.FDR.probes <- de_probes$probeID[idx.probe]
layout_c <- layoutConfig(
  vals = c(c(0,3), c(0,5), c(2,1), c(0,4)),
  nrow = 4,
  ncol = 2,
  lwid = c(0.18, 3.0),
  lhei = c(0.35, 0.03, 3, 0.6)
)
pdf('./figures/3c-heatmap_Ctrl_vs_FDR.pdf',
    width = 6,
    height = 10)
generateHeatMapPair(
  normal_probe_data = normal_probe_data,
  probe_ID = Ctrl.FDR.probes,
  idx = c(idx.Ctrl,idx.FDR),
  layout_c = layout_c
)
dev.off()
