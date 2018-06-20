rm(list=ls())
load('./data/normal_probe_data.RData')
de_probes <- read.csv(
  './data/DE_probes_annotated_pvals_fold_change.csv',
  stringsAsFactors = FALSE)
idx.probe <- de_probes$celiac_vs_FDR.pval_adjusted <= 0.05 & abs(de_probes$celiac_vs_FDR.fc) >= 1
FDR.CeD.probes <- de_probes$probeID[idx.probe]
idx.probe <- de_probes$celiac_vs_control.pval_adjusted <= 0.05 & abs(de_probes$celiac_vs_control.fc) >= 1
Ctrl.CeD.probes <- de_probes$probeID[idx.probe]
idx.probe <- de_probes$FDR_vs_control.pval_adjusted <= 0.05 & abs(de_probes$FDR_vs_control.fc) >= 1
Ctrl.FDR.probes <- de_probes$probeID[idx.probe]
de.probes <- unique(c(
  FDR.CeD.probes,
  Ctrl.CeD.probes,
  Ctrl.FDR.probes)
)
expr.data <- normal_probe_data$df[de.probes,]

#---
# perform PCA
expr.mat <- as.matrix(expr.data)
pca.expr <- prcomp(expr.mat)
summary.stats <- summary(pca.expr)
diagnosis <- normal_probe_data$targets
diagnosis[diagnosis == "celiac"] <- 'CeD'
diagnosis[diagnosis == "celiac"] <- 'CeD'
df.plot <- data.frame(PC1=pca.expr$rotation[,1],
                      PC2=pca.expr$rotation[,2],
                      Diagnosis=diagnosis)
library(ggplot2)
pdf('./figures/4-PCA_differentially-expressed-probes-among-diagnosis-groups.pdf',
width = 5, height = 5.3)
ggplot(df.plot,aes(x=PC1,y=PC2)) +
  geom_point(aes(col=Diagnosis,shape=Diagnosis),size=3,alpha=0.8) +
  xlab(paste("PC1 (",round(summary.stats$importance[2,1]*100,2),"%)",sep="")) +
  ylab(paste("PC2 (",round(summary.stats$importance[2,2]*100,2),"%)",sep="")) +
  scale_color_manual(
    values = c(
      '#217844ff', #green, CeD
      '#0088aaff', #blue, control
      '#aa4400ff') # brown, FDR
  ) +
  scale_shape_manual(
    values = c(
      0, # CeD
      1, # control
      2 # FDR
    )
  )
dev.off()
#-----
train.dat <- data.frame(diagnosis=diagnosis,t(expr.mat))
save(train.dat, file='./data/train.dat.RData')