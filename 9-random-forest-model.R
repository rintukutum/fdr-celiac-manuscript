rm(list=ls())
source('./func-utils.R')
load('./data/train.dat.RData')
#----
# following code will take time
set.seed(87569)
bt.var <- smartBoruta(
  x = train.dat[,-1],
  y = train.dat[,1],
  pValue = 0.001,
  maxRuns = 500)
save(bt.var,
     file = './data/bt.var.RData')
#----
bt.var_confirmed <- unlist(
  lapply(bt.var,
         function(x){
           names(x[x == 'Confirmed'])
           }
        )
  )
library(randomForest)
set.seed(10475)
rf <- randomForest(
  x=train.dat[,bt.var_confirmed],
  y=train.dat[,1],
  ntree = 10000,
  importance = TRUE,
  proximity = TRUE
)
write.csv(rf$importance,'./data/rf_variable_importance.csv')
rm(list=ls())
var_imp <- read.csv(
  './data/rf_variable_importance.csv',
  stringsAsFactors = FALSE,
  row.names = 1)
load('./data/train.dat.RData')
# build models using top 100
rf.imp <- var_imp[order(var_imp[,'MeanDecreaseAccuracy'],decreasing = TRUE),]
set.seed(1928)
library(randomForest)
rf.100 <- randomForest(
  x=train.dat[,rownames(rf.imp)[1:100]],
  y=train.dat[,1],
  ntree = 10000,
  importance = TRUE,
  proximity = TRUE
)
mds <- MDSplot(rf.100,train.dat[,1])
df.plot<- data.frame(mds$points,diagnosis =train.dat[,1])
library('ggplot2')
pdf('./figures/9-rf_100_top_imp_MDSplot.pdf',
  width = 5,height = 5.5)
ggplot(df.plot,aes(Dim.1,Dim.2)) +
  geom_point(aes(col=diagnosis,shape=diagnosis),size=3,alpha=0.8) +
  scale_color_manual(
    values = c(
      '#217844ff', #green
      '#0088aaff', #blue
      '#aa4400ff') # brown
    ) +
  xlab('Dim 1') +
  ylab('Dim 2') +
  scale_shape_manual(
    values = c(
      0,
      1,
      2
    )
  ) +
  ggtitle("MDS plot based on Random forests algo.\nwith top 100 features")
dev.off()
