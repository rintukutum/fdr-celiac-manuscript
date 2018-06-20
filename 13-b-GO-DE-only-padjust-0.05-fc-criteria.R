rm(list=ls())
de_probes <- read.csv(
  './manuscript/DE_probes_annotated_pvals_fold_change--only_padjust_0.05.csv',
  stringsAsFactors = FALSE
)
idx.fc.passed <- abs(de_probes$log2.FC) >=1

de_probes.comp <-plyr::dlply(
  de_probes[idx.fc.passed,],
  'comparison'
)
source('./func-utils.R')
go.comp <- list()
for(i in 1:length(de_probes.comp)){
  print(paste('i = ',i,sep=''))
  cond <- c('up','down')
  go.ind <- list()
  for(j in 1:length(cond)){
    print(j)
    entrezID.go <- getEntrezSig(
      de_probes.comp[[i]],
      type = cond[j])
    go <- doEnrichGO(x = entrezID.go)
    go.df <- plyr::ldply(go)
    colnames(go.df)[1] <- 'GO'
    go.ind[[j]] <- go.df
  }
  names(go.ind) <- cond
  go.ind.df <- plyr::ldply(go.ind)
  colnames(go.ind.df)[1] <- 'status'
  go.comp[[i]] <- go.ind.df
}
names(go.comp) <- names(de_probes.comp)
go.comp_fc_padjust <- go.comp
save(
  go.comp_fc_padjust,
  file = './data/go.comp_fc_padjust.RData'
)
go.comp.df <- plyr::ldply(
  go.comp_fc_padjust
)
colnames(go.comp.df)[1] <- 'Comparison'
write.csv(
  go.comp.df,
  file = './data/GO_DE_padjust_0.05_fc_1.csv',
  row.names = FALSE
)
