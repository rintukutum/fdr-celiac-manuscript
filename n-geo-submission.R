rm(list=ls())
load('data/preprocessing/celiac_data.RData');
load('data/preprocessing/celiac_control_data.RData');
load('data/preprocessing/celiac_family_data.RData')
#---
regular.probes <- function(x){
  idx <- which(x$genes$Status == 'regular')
  return(x[idx,])
}
cel <- regular.probes(x = celiac_data$raw)
con <- regular.probes(x = celiac_control_data$raw)
fdr <- regular.probes(x = celiac_family_data$raw)
common.raw.probes <- intersect(
  intersect(
    rownames(cel$E),
    rownames(con$E)
  ),
  rownames(fdr$E)
)

getRaw4GEO <- function(x,s.code){
  raw <- matrix(
    data = NA,
    nrow = length(common.raw.probes),
    ncol = 2*12)
  pos <- seq(1,24,by = 2)
  j <- 1
  for(i in pos){
    raw[,i] <- x$E[common.raw.probes,j]
    raw[,i+1] <- x$other$Detection[common.raw.probes,j]
    j <- j + 1
  }
  raw <- data.frame(raw)
  
  colnames(raw)[pos] <- s.code[colnames(x$E)]
  colnames(raw)[pos+1] <- 'Detection Pval'
  rownames(raw) <- common.raw.probes
  return(raw)
}
load('./data/sample_info.RData')
s.code <- gsub('[[:space:]]','',sample_info$Code)
names(s.code) <- sample_info$ID
cel.raw.geo <- getRaw4GEO(x = cel,
                          s.code)
con.raw.geo <- getRaw4GEO(x = con,
                          s.code)
fdr.raw.geo <- getRaw4GEO(x = fdr,
                          s.code)
raw.geo <- cbind(
  cel.raw.geo,
  con.raw.geo,
  fdr.raw.geo
)
write.csv(raw.geo,file = './geo-submission/raw-geo.csv')
df.scode <- data.frame(order = colnames(raw.geo)[grep('^C',colnames(raw.geo))],s.code)
write.csv(df.scode,
          file = './geo-submission/df_scode.csv')
##----
cel.norm <- regular.probes(x = celiac_data$normalized_quantile)
con.norm <- regular.probes(x = celiac_control_data$normalized_quantile)
fdr.norm <- regular.probes(x = celiac_family_data$normalized_quantile)

getNorm4GEO <- function(x,s.code){
  raw <- matrix(
    data = NA,
    nrow = length(common.raw.probes),
    ncol = 2*12)
  pos <- seq(1,24,by = 2)
  j <- 1
  for(i in pos){
    raw[,i] <- x$E[common.raw.probes,j]
    raw[,i+1] <- x$other$Detection[common.raw.probes,j]
    j <- j + 1
  }
  raw <- data.frame(raw)
  
  colnames(raw)[pos] <- s.code[colnames(x$E)]
  colnames(raw)[pos+1] <- 'Detection Pval'
  rownames(raw) <- common.raw.probes
  return(raw)
}

cel.norm.geo <- getRaw4GEO(x = cel.norm,
                           s.code)
con.norm.geo <- getRaw4GEO(x = con.norm,
                           s.code)
fdr.norm.geo <- getRaw4GEO(x = fdr.norm,
                           s.code)
norm.geo <- cbind(
  cel.norm.geo,
  con.norm.geo,
  fdr.norm.geo
)
write.csv(norm.geo,
          file = './geo-submission/norm-geo.csv')
####
#### for only common probes passed the criteria
#### filtered data
load('./data/common_probe_data.RData')
common.probes <- rownames(common_probe_data[[1]]$E)
filtered.data <- norm.geo[common.probes,]
write.csv(filtered.data,
          file = './geo-submission/filtered-data-geo.csv')

