rm(list=ls())
#-----
# load significant probe list among groups
de_probes <- read.csv(
  './manuscript/DE_probes_annotated_pvals_fold_change--only_padjust_0.05.csv',
  stringsAsFactors = FALSE
)
# fold change +/- log(FC) and padjust 0.05
idx.sig <- abs(de_probes$log2.FC) >= 1
de_probes.sig <- de_probes[idx.sig,]
probes <- unique(de_probes.sig$probeID)
#-----
# gene expression data 
load('./data/common_probe_data.RData')
de_probes.sig_data <- normal_probe_data
avg.expr.pheno <- lapply(common_probe_data,
                         function(x){
                           df <- data.frame(apply(x$E[probes,],1,mean))
                           colnames(df) <- x$name
                           return(df)
                         })
avg.expr.df <- do.call(cbind,
                       avg.expr.pheno)
pattern.expr <- apply(avg.expr.df,
                      1,
                      function(x){
                        y <- x
                        names(y) <- colnames(avg.expr.df)
                        y.out <- paste(names(sort(y)),collapse =' < ')
                        return(y.out)
                      })
load('./data/anno.de.df.RData')
rownames(anno.de.df) <- anno.de.df$probeID
pattern.df <- data.frame(
  probeID  = names(pattern.expr),
  EntrezID = anno.de.df[names(pattern.expr),'uid'],
  geneSymbol = anno.de.df[names(pattern.expr),'name'],
  pattern = pattern.expr,
  stringsAsFactors = FALSE
)
write.csv(pattern.df,
          './data/pattern-gene_expression-de-fc-criteria-1-padjust-0.05.csv',
          row.names = FALSE)
#---------------
rm(list=ls())
pattern.df <- read.csv(
  './data/pattern-gene_expression-de-fc-criteria-1-padjust-0.05.csv',
  stringsAsFactors = FALSE
)
pattern.df.group <- plyr::dlply(
  pattern.df,
  'pattern'
)
pattern.df.group.genes <- sapply(
  pattern.df.group,
  function(x){
    length(unique(na.omit(x$EntrezID)))
  }
)
pdf('./figures/14-pattern-analysis-stats.pdf',
    width = 4,
    height = 5)
barplot(
  pattern.df.group.genes,
  axisnames = FALSE,
  ylab = '# Annotated genes',
  col = c(
    '#806600ff',
    '#d4aa00ff',
    '#decd87ff',
    '#c8c4b7ff'),
  las = 2
)
legend(
  x = 0.2,
  y = 700,
  legend = c(
      expression(mu['(CeD)']~'<'~mu['(Control)']~'<'~mu['(FDR)']),
      expression(mu['(CeD)']~'<'~mu['(FDR)']~'<'~mu['(Control)']),
      expression(mu['(FDR)']~'<'~mu['(CeD)']~'<'~mu['(Control)']),
      expression(mu['(FDR)']~'<'~mu['(Control)']~'<'~mu['(CeD)'])
  ),
  fill =  c(
    '#806600ff',
    '#d4aa00ff',
    '#decd87ff',
    '#c8c4b7ff'),
  cex = 0.85
)
dev.off()
###############
rm(list=ls())
pattern.df <- read.csv(
  './data/pattern-gene_expression-de-fc-criteria-1-padjust-0.05.csv'
)
pattern.df.group <- plyr::dlply(
  pattern.df,
  'pattern'
)
source('./func-utils.R')
go.pattern <- list()
for(i in 1:length(pattern.df.group)){
  print(i)
  go.pattrn <- doEnrichGO(
    x = na.omit(pattern.df.group[[i]]$EntrezID)
  )
  go.pattrn.df <- plyr::ldply(go.pattrn)
  colnames(go.pattrn.df)[1] <- 'GO'
  go.pattern[[i]] <- go.pattrn.df
}
names(go.pattern) <- names(pattern.df.group)
go.pattern.df <- plyr::ldply(go.pattern)
colnames(go.pattern.df)[1] <- 'pattern'
go.pattern.df.cutoff <- go.pattern.df
save(
  go.pattern.df.cutoff,
  file = './data/go.pattern.df.cutoff.RData'
)
write.csv(
  go.pattern.df.cutoff,
  file = './data/GO_pattern-cutoff-adjusted-0.05_de_padjust_fc_criteria.csv',
  row.names = FALSE
)
##########################
go.pattern <- list()
for(i in 1:length(pattern.df.group)){
  print(i)
  go.pattrn <- doEnrichGO(
    x = na.omit(pattern.df.group[[i]]$EntrezID),
    pval.cutoff = 1,
    qval.cutoff = 1
  )
  go.pattrn.df <- plyr::ldply(go.pattrn)
  colnames(go.pattrn.df)[1] <- 'GO'
  go.pattern[[i]] <- go.pattrn.df
}
names(go.pattern) <- names(pattern.df.group)
go.pattern.df <- plyr::ldply(go.pattern)
colnames(go.pattern.df)[1] <- 'pattern'
############
go.pattern.df.no.cutoff <- go.pattern.df
idx.sig <- go.pattern.df.no.cutoff$pvalue <=0.05
go.pattern.df.no.cutoff <- go.pattern.df.no.cutoff[idx.sig,]
save(
  go.pattern.df.no.cutoff,
  file = './data/go.pattern.df.no.cutoff.RData'
)
write.csv(
  go.pattern.df.no.cutoff,
  file = './data/GO_pattern-cutoff-pval-0.05_de_padjust_fc_criteria.csv',
  row.names = FALSE
)
#############
rm(list=ls())
#### 
#### show top 10 processes
load('./data/go.pattern.df.no.cutoff.RData')
pattern.go <- plyr::dlply(
  go.pattern.df.no.cutoff,
  .variables = 'GO',
  .fun = function(x){
    x.l <- plyr::dlply(
      x,
      .variables = 'pattern',
      .fun = function(x){
        x[order(x$pvalue,decreasing = FALSE)[1:20],]
      })
    return(x.l)
  }
)
library(ggplot2)

exprs.pattern <- c(
  'celiac < control < FDR' = expression(mu['(CeD)']~'<'~mu['(Control)']~'<'~mu['(FDR)']),
  'celiac < FDR < control' = expression(mu['(CeD)']~'<'~mu['(FDR)']~'<'~mu['(Control)']),
  'FDR < celiac < control' = expression(mu['(FDR)']~'<'~mu['(CeD)']~'<'~mu['(Control)']),
  'FDR < control < celiac' = expression(mu['(FDR)']~'<'~mu['(Control)']~'<'~mu['(CeD)'])
)
p <- list()
jj <- 1
dir.create('./figures/pattern/',showWarnings = FALSE)
for(i in 1:length(pattern.go)){
  for(j in 1:length(pattern.go[[i]])){
    dat <- pattern.go[[i]][[j]]
    pattern.name <- exprs.pattern[unique(dat$pattern)]
    p[[jj]] <- ggplot(dat,
                      aes(y=Description,x=pattern)) +
      geom_tile(aes(fill=-log10(pvalue)),
                width=0.95,
                height=0.9) +
      xlab('') +
      scale_fill_gradient(limits=c(1,10),low = '#ffd5d5ff',high = '#ff2a2aff',guide=guide_colorbar(raster=FALSE)) +
      scale_y_discrete(limits=rev(dat$Description)) +
      scale_x_discrete(labels = pattern.name) +
      facet_wrap(facets = 'GO')
    
    out.name <- paste(c(unique(dat$pattern),unique(dat$GO)),collapse = '--')
    pdf(paste('./figures/pattern/go--',out.name,'.pdf',sep=''),
        width = 6.5,
        height = 6,
        useDingbats = FALSE)
    print(p[[jj]])
    dev.off()
    jj <- jj + 1
  }
}
