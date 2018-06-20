rm(list=ls());
dat.anno <- read.delim(
  file='./data/annotation/HumanHT-12_V4_0_R1_15002873_B_extracted.txt',
  header=TRUE,
  stringsAsFactors = FALSE,
  sep = '\t',
  row.names = NULL)

probe2genes_data<-dat.anno[,c("Probe_Id",
                              "ILMN_Gene",
                              "Entrez_Gene_ID",
                              "Symbol")];
idx.no.symbol <- sapply(probe2genes_data$Symbol,nchar) == 0
probe2genes_data$Symbol[idx.no.symbol] <- probe2genes_data$ILMN_Gene[idx.no.symbol]
save(probe2genes_data,file='./data/probe2genes_data.RData')
rm(list=ls());
load('./data/fc_pval_pair_wise_test.RData')
load('./data/probe2genes_data.RData')
gene.entrez <- data.frame(
  probeID = probe2genes_data$Probe_Id,
  entrezID = probe2genes_data$Entrez_Gene_ID,
  stringsAsFactors = FALSE
)

gene.de <- gene.entrez[gene.entrez$probeID %in% rownames(fc_pval_pair_wise_test),]
source('./func-utils.R')
anno.de <- getDetailGeneInfo(
  entrezIDs = gene.de$entrezID
)
library(plyr)
anno.de.df <- ddply(
  .data = gene.de,
  .variables = 'probeID',
  .fun = function(x){
    anno.de[as.character(x$entrezID),1:5]
  }
)
for(i in 1:nrow(anno.de.df)){
  if(!is.na(anno.de.df$currentid[i])){
  if(nchar(anno.de.df$currentid[i]) == 0){
    anno.de.df$currentid[i] <- anno.de.df$uid[i]
  }
  }
}
save(anno.de.df,
     file = './data/anno.de.df.RData')
#-------------
load('./data/anno.de.df.RData')
load('./data/probe2genes_data.RData')
idx.no.symbol <- is.na(anno.de.df$name)
probes <- anno.de.df$probeID[idx.no.symbol]
rownames(anno.de.df) <- anno.de.df$probeID

rownames(probe2genes_data) <- probe2genes_data$Probe_Id
anno.de.df[probes,'name'] <- probe2genes_data[probes,'Symbol']
write.csv(
  anno.de.df,
  './data/differentially_expressed_probes_annotation_entrez.csv',
  row.names = FALSE
)
anno.de.df <- anno.de.df[rownames(fc_pval_pair_wise_test),]
de.probes.pval.fc.anno <- cbind(fc_pval_pair_wise_test,
                                anno.de.df[,c('probeID','currentid','name','description')])
de.probes.pval.fc.anno <- de.probes.pval.fc.anno[,c(10:12,1:9,13)]
write.csv(
  de.probes.pval.fc.anno,
  './data/DE_probes_annotated_pvals_fold_change.csv',
  row.names = FALSE
)
####
rm(list=ls())
de.probes <- read.csv(
  './data/DE_probes_annotated_pvals_fold_change.csv',
  stringsAsFactors = FALSE
)
de.ced_fdr <- de.probes[,c(1:3,13,grep('celiac_vs_FDR',colnames(de.probes)))]

colnames(de.ced_fdr)[5:7] <- c('p.value','log2.FC','p.adjust')
idx.sig <- de.ced_fdr$p.adjust <= 0.05
de.ced_fdr <- cbind(
  de.ced_fdr[idx.sig,],
  comparison = 'CeD/FDR'
)
de.ced_ctrl <- de.probes[,c(1:3,13,grep('celiac_vs_control',colnames(de.probes)))]
colnames(de.ced_ctrl)[5:7] <- c('p.value','log2.FC','p.adjust')
idx.sig <- de.ced_ctrl$p.adjust <= 0.05
de.ced_ctrl <- cbind(
  de.ced_ctrl[idx.sig,],
  comparison = 'CeD/Control'
)
de.fdr_ctrl <- de.probes[,c(1:3,13,grep('FDR_vs_control',colnames(de.probes)))]
colnames(de.fdr_ctrl)[5:7] <- c('p.value','log2.FC','p.adjust')
idx.sig <- de.fdr_ctrl$p.adjust <= 0.05
de.fdr_ctrl <- cbind(
  de.fdr_ctrl[idx.sig,],
  comparison = 'FDR/Control'
)
de.all <- rbind(
  de.ced_fdr,
  de.ced_ctrl,
  de.fdr_ctrl
)
colnames(de.all)[2] <- 'EntrezID'
write.csv(
  de.all[,-grep('description', colnames(de.all))],
  './manuscript/DE_probes_annotated_pvals_fold_change--only_padjust_0.05.csv',
  row.names = FALSE
)
de.probes.details <- de.all[!duplicated(de.all$probeID),1:4]
write.csv(
  de.probes.details,
  './manuscript/DE_probes_annotation--only_padjust_0.05.csv',
  row.names = FALSE
)
